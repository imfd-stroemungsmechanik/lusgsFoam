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
    U_.ref() =
        rhoU_()
       /rho_();
    U_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() == rho_.boundaryField()*U_.boundaryField();

    volScalarField muEff("muEff", turbulence_.muEff());

    // Solve momentum equation
    if (!inviscid_)
    {
        solve
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
          - fvm::laplacian(muEff, U_)
          - fvc::div(flux_.tauMC())
        );
        rhoU_ = rho_*U_;
    }

    // Update internal energy
    thermo_.he() = rhoE_/rho_ - 0.5*magSqr(U_);
    thermo_.he().correctBoundaryConditions();
    thermo_.correct();
    rhoE_.boundaryFieldRef() ==
        rho_.boundaryField()*
        (
            thermo_.he().boundaryField() + 0.5*magSqr(U_.boundaryField())
        );

    // Solve energy equation
    if (!inviscid_)
    {
        solve
        (
            fvm::ddt(rho_, thermo_.he()) - fvc::ddt(rho_, thermo_.he())
          - fvm::laplacian(turbulence_.alphaEff(), thermo_.he())
        );
        thermo_.correct();
        rhoE_ = rho_*(thermo_.he() + 0.5*magSqr(U_));
    }

    // Update pressure field
    p_.ref() =
        rho_()
       /thermo_.psi();
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
    omega_(1.0),
    inviscid_(true),
    gamma_(linearInterpolate(thermo_.Cp()/thermo_.Cv())),
    D_(mesh_.nCells(), Zero),
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
    deltaWStarRho_
    (
        IOobject
        (
            "lusgs::deltaWStarRho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaWStarRho", dimDensity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWStarRhoU_
    (
        IOobject
        (
            "lusgs::deltaWStarRhoU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("deltaWStarRho", dimDensity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWStarRhoE_
    (
        IOobject
        (
            "lusgs::deltaWStarRhoE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaWStarRho", dimDensity*dimVelocity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWRho_
    (
        IOobject
        (
            "lusgs::deltaWRho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_
    ),
    deltaWRhoU_
    (
        IOobject
        (
            "lusgs::deltaWRhoU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoU_
    ),
    deltaWRhoE_
    (
        IOobject
        (
            "lusgs::deltaWRhoE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoE_
    )
{
    Info<< "\nInitializing LU-SGS scheme" << endl
        << "    omega = " << omega_ << endl;

    if (max(thermo_.mu()().primitiveField()) > 0.0)
    {
        inviscid_ = false;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::LUSGS::sweep()
{
    // Add time derivative to residuals
    const volScalarField R_rho
    (
        fvc::ddt(rho_)
      + fvc::div(flux_.phi())
    );

    const volVectorField R_rhoU
    (
        fvc::ddt(rhoU_)
      + fvc::div(flux_.phiUp())
      - fvc::div(flux_.tauMC())
    );
    const volScalarField R_rhoE
    (
        fvc::ddt(rhoE_)
      + fvc::div(flux_.phiEp())
      - fvc::div(flux_.sigmaDotU())
    );

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
        dimensionedScalar("lusgs::meshPhi", dimVolume/dimTime, Zero)
    );
    if (mesh_.moving())
    {
        meshPhi = mesh_.phi();
    }

    // Time step
    const scalar dt = mesh_.time().deltaTValue();

    // Fill diagonal matrix D with cell values
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

    // Loop over all patches
    forAll(mesh_.boundary(), patchI)
    {
        const labelUList& pFaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& pSf = mesh_.Sf().boundaryField()[patchI];
        const scalarField& pMagSf = mesh_.magSf().boundaryField()[patchI];
        const scalarField& pMeshPhi = meshPhi.boundaryField()[patchI];
        const vectorField& pCf = mesh_.Cf().boundaryField()[patchI];

        // Loop over all faces on patchI
        forAll(mesh_.boundary()[patchI], faceI)
        {
            label cellI = pFaceCells[faceI];

            //scalar dvol = mag( delta[faceI] & pSf[faceI] );
            scalar dvol = mag( (mesh_.C()[cellI] - pCf[faceI]) & pSf[faceI] );

            // Calculate spectral radii
            scalar ac = 0.5 * omega_ * ( mag((U_[cellI] & pSf[faceI]) - pMeshPhi[faceI])
              + a[cellI]*mag(pSf[faceI]));
            scalar av = sqr(pMagSf[faceI]) / dvol * nuMax[cellI];

            D_[cellI] += ac + av;
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
        scalar L_rho = Zero;
        vector L_rhoU = Zero;
        scalar L_rhoE = Zero;

        // Loop over all faces of cellI
        forAll(cellFaces, i)
        {
            // Face ID
            const label faceI = cellFaces[i];

            // Get owner and neighbour cell ID
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            // Check if face is internal
            if (mesh_.isInternalFace(faceI))
            {
                const scalar dvol = mag( (mesh_.C()[own] - mesh_.C()[nei]) & Sf[faceI] );

                // Calculate spectral radii
                const scalar ac = 0.5 * omega_ * ( mag((U_[cellI] & Sf[faceI]) - meshPhi[faceI] ) + a[cellI] * magSf[faceI] );
                const scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];
		
                // Calculate lower-diagonal matrix L
                if (cellI == nei)
                {
                    const scalar rho1 = rho_[own] + deltaWStarRho_[own];
                    const vector rhoU1 = rhoU_[own] + deltaWStarRhoU_[own];
                    const scalar rhoE1 = rhoE_[own] + deltaWStarRhoE_[own];

                    const scalar p1 = (gamma_[faceI] - 1.0) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    // Previous and new fluxes
                    const scalar phi0 = U_[own] & Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1) & Sf[faceI];

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
        deltaWStarRho_[cellI] = (-R_rho[cellI]*V[cellI] + L_rho) / D_[cellI];
        deltaWStarRhoU_[cellI] = (-R_rhoU[cellI]*V[cellI] + L_rhoU) / D_[cellI];
        deltaWStarRhoE_[cellI] = (-R_rhoE[cellI]*V[cellI] + L_rhoE) / D_[cellI];
    }

    //
    // Backward sweep
    //

    deltaWRho_ *= 0.0;
    deltaWRhoU_ *= 0.0;
    deltaWRhoE_ *= 0.0;

    // Reverse loop over all cells
    forAllReverse(mesh_.cells(), cellI)
    {
        // List of faces of cellI
        const labelList& cellFaces = mesh_.cells()[cellI];

        // Sum of upper-diagonal matrix
        scalar U_rho = Zero;
        vector U_rhoU = Zero;
        scalar U_rhoE = Zero;

        // Loop over all faces of cellI
        forAll(cellFaces, i)
        {
            // Face ID
            const label faceI = cellFaces[i];

            // Get owner and neighbour cell ID
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            // Check if face is internal
            if (mesh_.isInternalFace(faceI))
            {
                // Calculate upper-diagonal matrix U
                if (cellI == own)
                {
                    const scalar dvol = mag( (mesh_.C()[own] - mesh_.C()[nei]) & Sf[faceI] );

                    // Calculate spectral radii
                    const scalar ac = 0.5 * omega_ * ( mag((U_[cellI] & Sf[faceI]) - meshPhi[faceI]) + a[cellI] * magSf[faceI] );
                    const scalar av = sqr(magSf[faceI]) / dvol * nuMax[nei];

                    const scalar rho1 = rho_[nei] + deltaWRho_[nei];
                    const vector rhoU1 = rhoU_[nei] + deltaWRhoU_[nei];
                    const scalar rhoE1 = rhoE_[nei] + deltaWRhoE_[nei];

                    const scalar p1 = (gamma_[faceI] - 1.0) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    // Previous and new fluxes
                    const scalar phi0 = U_[nei] & Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1) & Sf[faceI];

                    // Relative fluxes
                    const scalar phi0r = phi0 - meshPhi[faceI];
                    const scalar phi1r = phi1 - meshPhi[faceI];

                    // Speed of sound at interface
                    U_rho += (ac + av)*deltaWRho_[nei]
                      - 0.5*(rho1*phi1r - rho_[nei]*phi0r);
                    U_rhoU += (ac + av)*deltaWRhoU_[nei]
                      - 0.5*(rhoU1*phi1r - rhoU_[nei]*phi0r + (p1 - p_[nei])*Sf[faceI]);
                    U_rhoE += (ac + av)*deltaWRhoE_[nei]
                      - 0.5*((rhoE1*phi1r + p1*phi1) - (rhoE_[nei]*phi0r + p_[nei]*phi0));
                }
            }
        }

        // Backwards sweep
        deltaWRho_[cellI] = deltaWStarRho_[cellI] + U_rho / D_[cellI];
        deltaWRhoU_[cellI] = deltaWStarRhoU_[cellI] + U_rhoU / D_[cellI];
        deltaWRhoE_[cellI] = deltaWStarRhoE_[cellI] + U_rhoE / D_[cellI];

    }

    // Update conserative fields
    rho_ += deltaWRho_;
    rhoU_ += deltaWRhoU_;
    rhoE_ += deltaWRhoE_;

    // Update primitive fields
    updatePrimitiveFields();

}

// ************************************************************************* //
