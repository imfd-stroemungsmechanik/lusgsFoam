/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA


\*---------------------------------------------------------------------------*/

#include "LUSGS.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::LUSGS::updatePrimitiveFields()
{
    // Update velocity
    U_.primitiveFieldRef() =
        rhoU_.primitiveField()
       /rho_.primitiveField();
    U_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() == rho_.boundaryField()*U_.boundaryField();

    // Update internal energy
    thermo_.he() = rhoE_/rho_ - 0.5*magSqr(U_);
    thermo_.he().correctBoundaryConditions();
    thermo_.correct();
    rhoE_.boundaryFieldRef() ==
        rho_.boundaryField()*
        (
            thermo_.he().boundaryField() + 0.5*magSqr(U_.boundaryField())
        );

    // Update pressure field
    p_.primitiveFieldRef() =
        rho_.primitiveField()
       /thermo_.psi().primitiveField();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() == thermo_.psi().boundaryField()*p_.boundaryField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::LUSGS::LUSGS
(
    volScalarField& p,
    volScalarField& rho,
    volVectorField& U,
    numericFlux& flux,
    psiThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    mesh_(rho.mesh()),
    p_(p),
    rho_(rho),
    U_(U),
    flux_(flux),
    thermo_(thermo),
    turbulence_(turbulence),
    omega_(
        max
        (
            1.0,
            min
            (
                2.0,
                mesh_.solutionDict().subDict("LUSGS").lookupOrDefault<scalar>("omega", scalar(1.0))
            )
        )
    ),
    maxIter_(readLabel(mesh_.solutionDict().subDict("LUSGS").lookup("maxIter"))),
    relTol_(readScalar(mesh_.solutionDict().subDict("LUSGS").lookup("relTol"))),
    gamma_(linearInterpolate(thermo_.Cp()/thermo_.Cv())),
    D_(mesh_.nCells(), scalar(0.0)),
    rhoU_
    (
        IOobject
        (
            "lusgs::rhoU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ),
    rhoE_
    (
        IOobject
        (
            "lusgs::rhoE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*(thermo_.he() + 0.5*magSqr(U_))
    ),
    deltaWStarRho_(rho_),
    deltaWStarRhoU_(rhoU_),
    deltaWStarRhoE_(rhoE_),
    resRho_(0.0),
    resRhoU_(0.0),
    resRhoE_(0.0)
{
    Info<< "\nInitializing LU-SGS scheme" << endl
        << "    omega = " << omega_ << endl
        << "    maxIter = " << maxIter_ << endl
        << "    relTol = " << relTol_ << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::LUSGS::solve()
{
    // Compute inviscid fluxes
    flux_.update();

    // Add time derivative to residuals
    const volScalarField R_rho = flux_.resRho() + fvc::ddt(rho_);
    const volVectorField R_rhoU = flux_.resRhoU() + fvc::ddt(rhoU_);
    const volScalarField R_rhoE = flux_.resRhoE() + fvc::ddt(rhoE_);

    // Owner and neighbour face index
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    // Face area vector and cell volume
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const scalarField& V = mesh_.V();

    // Flux of mesh
    surfaceScalarField meshPhi
    (
        IOobject
        (
            "lusgs::meshPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("lusgs::meshPhi", dimVolume/dimTime, scalar(0.0))
    );
    if (mesh_.moving())
    {
        meshPhi = mesh_.phi();
    }

    // Time step
    const scalar dt = mesh_.time().deltaTValue();

    // Fill diagonal matrix Dref with cell values
    forAll(mesh_.cells(), cellI)
    {
        // Calculate first part of D
        D_[cellI] = V[cellI] / dt;
    }

    // Speed of sound
    volScalarField a = sqrt(thermo_.Cp()/thermo_.Cv() / thermo_.psi());

    // Viscous spectral radius
    volScalarField nuMax
    (
        max
        (
            4.0/3.0*turbulence_.nuEff(),
            thermo_.Cp()/thermo_.Cv()*turbulence_.alphaEff()/rho_
        )
    );

    // Loop over all boundaries
    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh_.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh_.Sf().boundaryField()[patchi];
        const scalarField& pMeshPhi = meshPhi.boundaryField()[patchi];

        forAll(mesh_.boundary()[patchi], faceI)
        {

            label cellI = pFaceCells[faceI];
            scalar ac = 0.5 * omega_ * ( mag((U_[cellI] & pSf[faceI]) - pMeshPhi[faceI])
              + a[cellI]*mag(pSf[faceI]));

            D_[cellI] += ac;
        }
    }

    // Reset the values for deltaWStar
    deltaWStarRho_ *= 0.0;
    deltaWStarRhoU_ *= 0.0;
    deltaWStarRhoE_ *=0.0;



    //
    // Forward sweep
    //

    // Forward loop over all cells
    forAll(mesh_.cells(), cellI)
    {
        // List of faces of cellI
        const labelList& cellFaces = mesh_.cells()[cellI];

        // Sum of lower-diagonal matrix
        scalar L_rho = 0.0;
        vector L_rhoU = vector(0,0,0);
        scalar L_rhoE = 0.0;

        // Loop over all faces of cellI
        forAll(cellFaces, i)
        {
            // Face ID
            label faceI = cellFaces[i];

            // Check if face is internal
            if (mesh_.isInternalFace(faceI))
            {
                // Get owner and neighbour cell ID
                label own = owner[faceI];
                label nei = neighbour[faceI];

                scalar dvol = mag( (mesh_.C()[own] - mesh_.C()[nei]) & Sf[faceI] );

                scalar ac = 0.5 * omega_ * ( mag((U_[cellI] & Sf[faceI]) - meshPhi[faceI] ) + a[cellI] * magSf[faceI] );
                scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];
		
                // Calculate L
                if (cellI == nei)
                {
                    scalar rho1 = rho_[own] + deltaWStarRho_[own];
                    vector rhoU1 = rhoU_[own] + deltaWStarRhoU_[own];
                    scalar rhoE1 = rhoE_[own] + deltaWStarRhoE_[own];

                    scalar p1 = (gamma_[faceI] - 1.0) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    // Previous and new fluxes
                    scalar phi0 = U_[own] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];

                    // Relative fluxes
                    const scalar phi0r = phi0 - meshPhi[faceI];
                    const scalar phi1r = phi1 - meshPhi[faceI];

                    // Calculate lower diagonal matrix
                    L_rho += (ac + av)*deltaWStarRho_[own]
                      + 0.5*(rho1*phi1r - rho_[own]*phi0r);
                    L_rhoU += (ac + av)*deltaWStarRhoU_[own]
                      + 0.5*(rhoU1*phi1r - rhoU_[own]*phi0r + (p1 - p_[own])*Sf[faceI]);
                    L_rhoE += (ac + av)*deltaWStarRhoE_[own]
                      + 0.5*((rhoE1*phi1r + p1*phi1) - (rhoE_[own]*phi0r + p_[own]*phi0));
                }

                D_[cellI] += ac + av;
            }
        }

        // Perform calculation for deltaWStar
        deltaWStarRho_[cellI]  = ( -R_rho[cellI]*V[cellI] +  L_rho) / D_[cellI];
        deltaWStarRhoU_[cellI] = (-R_rhoU[cellI]*V[cellI] + L_rhoU) / D_[cellI];
        deltaWStarRhoE_[cellI] = (-R_rhoE[cellI]*V[cellI] + L_rhoE) / D_[cellI];
    }



    //
    // Backward sweep
    //

    volScalarField deltaWRho(0.*deltaWStarRho_);
    volVectorField deltaWRhoU(0.*deltaWStarRhoU_);
    volScalarField deltaWRhoE(0.*deltaWStarRhoE_);

    // Reverse loop over all cells
    forAllReverse(mesh_.cells(), cellI)
    {
        // List of faces of cellI
        const labelList& cellFaces = mesh_.cells()[cellI];

        // Sum of upper-diagonal matrix
        scalar U_rho = 0.0;
        vector U_rhoU = vector(0,0,0);
        scalar U_rhoE = 0.0;

        // Loop over all faces of cellI
        forAll(cellFaces, i)
        {
            // Face ID
            label faceI = cellFaces[i];

            // Check if face is internal
            if (mesh_.isInternalFace(faceI))
            {
                // Get owner and neighbour cell ID
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // Calculate U
                if (cellI == own)
                {
                    scalar dvol = mag( (mesh_.C()[own] - mesh_.C()[nei]) & Sf[faceI] );
                    
                    scalar ac = 0.5 * omega_ * ( mag((U_[nei] & Sf[faceI]) - meshPhi[faceI]) + a[nei] * magSf[faceI] );
                    scalar av = sqr(magSf[faceI]) / dvol * nuMax[nei];

                    scalar rho1 = rho_[nei] + deltaWRho[nei];
                    vector rhoU1 = rhoU_[nei] + deltaWRhoU[nei];
                    scalar rhoE1 = rhoE_[nei] + deltaWRhoE[nei];

                    scalar p1 = (gamma_[faceI] - 1.0) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    // Previous and new fluxes
                    scalar phi0 = U_[nei] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];

                    // Relative fluxes
                    scalar phi0r = phi0 - meshPhi[faceI];
                    scalar phi1r = phi1 - meshPhi[faceI];

                    // Speed of sound at interface
                    U_rho += (ac + av)*deltaWRho[nei]
                      - 0.5*(rho1*phi1r - rho_[nei]*phi0r);
                    U_rhoU += (ac + av)*deltaWRhoU[nei]
                      - 0.5*(rhoU1*phi1r - rhoU_[nei]*phi0r + (p1 - p_[nei])*Sf[faceI]);
                    U_rhoE += (ac + av)*deltaWRhoE[nei]
                      - 0.5*((rhoE1*phi1r + p1*phi1) - (rhoE_[nei]*phi0r + p_[nei]*phi0));
                }
            }
        }

        // Backwards sweep
        deltaWRho[cellI] = deltaWStarRho_[cellI] + U_rho / D_[cellI];
        deltaWRhoU[cellI] = deltaWStarRhoU_[cellI] + U_rhoU / D_[cellI];
        deltaWRhoE[cellI] = deltaWStarRhoE_[cellI] + U_rhoE / D_[cellI];
    }

    // Update conserative fields
    rho_ += deltaWRho;
    rhoU_ += deltaWRhoU;
    rhoE_ += deltaWRhoE;

    // Calculate convergence criteria
    resRho_  = fvc::domainIntegrate( mag(deltaWRho) / dt ).value();
    resRhoU_ = fvc::domainIntegrate( mag(deltaWRhoU) / dt ).value();
    resRhoE_ = fvc::domainIntegrate( mag(deltaWRhoE) / dt ).value();

    // Update primitive fields
    updatePrimitiveFields();
}

// ************************************************************************* //
