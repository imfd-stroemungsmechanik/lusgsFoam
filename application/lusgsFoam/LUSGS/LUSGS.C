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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::LUSGS::updatePrimitiveFields()
{
    // Update velocity
    U_.ref() =
        rhoU_()
       /rho_();
    U_.correctBoundaryConditions();
    MRF_.correctBoundaryVelocity(U_);  
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
    volVectorField& U,
    volScalarField& rho,
    inviscidFlux& flux,
    psiThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const IOMRFZoneList& MRF
)
:
    mesh_(rho.mesh()),
    p_(p),
    U_(U),
    rho_(rho),
    flux_(flux),
    thermo_(thermo),
    turbulence_(turbulence),
    MRF_(MRF),
    gamma_(linearInterpolate(thermo_.gamma())),
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
            mesh_
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
            mesh_
        ),
        mesh_,
        dimensionedVector("deltaWStarRhoU", dimDensity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWStarRhoE_
    (
        IOobject
        (
            "lusgs::deltaWStarRhoE",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("deltaWStarRhoE", dimDensity*dimVelocity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWRho_
    (
        IOobject
        (
            "lusgs::deltaWRho",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("deltaWRho", dimDensity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWRhoU_
    (
        IOobject
        (
            "lusgs::deltaWRhoU",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector("deltaWRhoU", dimDensity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWRhoE_
    (
        IOobject
        (
            "lusgs::deltaWRhoE",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("deltaWRhoE", dimDensity*dimVelocity*dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    omega_(1.0),
    nIter_(1)
{
    Info<< "\nInitializing LU-SGS scheme" << endl;

    if (mesh_.solutionDict().isDict("LUSGS"))
    {
        omega_ = mesh_.solutionDict().subDict("LUSGS").lookupOrDefault<scalar>("omega", 1.0);
        nIter_ = mesh_.solutionDict().subDict("LUSGS").lookupOrDefault<label>("nIter", 1);
    }
    Info << "    omega       " << omega_ << endl;
    Info << "    nIter       " << nIter_ << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::LUSGS::sweep()
{
    // LUSGS iteration loop
    for (int iter = 0; iter < nIter_; iter++)
    {
        // Update flux after first iteration
        if (iter > 0)
        {
            flux_.update();
        }

        // Viscous and Reynolds stress tensor
        const volSymmTensorField tau
        (
            "tau",
          - turbulence_.devRhoReff()
          - ((2.0/3.0)*I)*rho_*turbulence_.k()
        );

        // Assemble residuals including time derivatives
        const volScalarField R_rho
        (
            fvc::ddt(rho_)
          + fvc::div(flux_.phi())
        );
        const volVectorField R_rhoU
        (
            fvc::ddt(rhoU_)
          + fvc::div(flux_.phiUp())
          - fvc::div(tau) 
          + MRF_.DDt(rho_, U_)
        );
        const volScalarField R_rhoE
        (
            fvc::ddt(rhoE_)
          + fvc::div(flux_.phiEp())
          - fvc::div(tau & U_, "div(tau&U)")
          - fvc::laplacian
            (
                (turbulence_.mu()+0.6*turbulence_.mut()),
                turbulence_.k()
            )
          - fvc::laplacian(turbulence_.alphaEff(), thermo_.he())
        );

        // Owner and neighbour face index
        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();

        // Face area vectors, face normal vector, and cell volume
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
        meshPhi.setOriented(true);

        // Compute mesh flux due to mesh motion ...
        if (mesh_.moving())
        {
            meshPhi = fvc::meshPhi(U_);
        }
        // ... or due to MRF
        MRF_.makeAbsolute(meshPhi);

        // Time step
        const scalar dt = mesh_.time().deltaTValue();

        // Fill diagonal matrix D with cell values
        forAll(mesh_.cells(), cellI)
        {
            // Calculate first part of D
            D_[cellI] = V[cellI] / dt;
        }

        // Speed of sound
        volScalarField c = sqrt(thermo_.gamma()*thermo_.p()/thermo_.rho());

        // Viscous spectral radius
        volScalarField nuMax
        (
            max
            (
                4.0/(3.0*rho_),
                thermo_.gamma() / rho_
            )
          * turbulence_.alphaEff()
        );

        // Initial residuals
        scalar initialRes_rho  = fvc::domainIntegrate( mag(R_rho) / dt ).value();
        scalar initialRes_rhoU = fvc::domainIntegrate( mag(R_rhoU) / dt ).value();
        scalar initialRes_rhoE = fvc::domainIntegrate( mag(R_rhoE) / dt ).value();

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

                scalar dvol = mag( (mesh_.C()[cellI] - pCf[faceI]) & pSf[faceI] );

                // Calculate spectral radii
                scalar lc = mag((U_[cellI] & pSf[faceI]) - pMeshPhi[faceI])
                  + c[cellI]*mag(pSf[faceI]);
                scalar lv = sqr(pMagSf[faceI]) / dvol * nuMax[cellI];

                D_[cellI] += 0.5*lc + lv;
            }
        }

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
    
                    // Arithmetic average
                    const vector Uf = 0.5*(U_[own] + U_[nei]);
                    const scalar cf = 0.5*(c[own] + c[nei]);
                    const scalar nuMaxf = 0.5*(nuMax[own] + nuMax[nei]);
    
                    // Calculate spectral radii
                    const scalar lc = mag( (Uf & Sf[faceI]) - meshPhi[faceI] )
                      + cf*magSf[faceI];
                    const scalar lv = sqr(magSf[faceI]) / dvol * nuMaxf;
    
                    D_[cellI] += 0.5*lc + lv;

                    // Calculate lower-diagonal matrix L
                    if (cellI == nei)
                    {
    
                        // Cell values of owner cell
                        const scalar rho = rho_[own] + deltaWStarRho_[own];
                        const vector rhoU = rhoU_[own] + deltaWStarRhoU_[own];
                        const scalar rhoE = rhoE_[own] + deltaWStarRhoE_[own];
    
                        const vector U = rhoU / rho;
                        const vector n = Sf[faceI] / magSf[faceI];
    
                        const scalar Ux = U.component(0);
                        const scalar Uy = U.component(1);
                        const scalar Uz = U.component(2);
                        const scalar nx = n.component(0);
                        const scalar ny = n.component(1);
                        const scalar nz = n.component(2);
    
                        // coefficients
                        const scalar phi = 0.5*(gamma_[faceI] - 1.0)*magSqr(U);
                        const scalar a1 = gamma_[faceI] * rhoE / rho - phi;
                        const scalar a2 = gamma_[faceI] - 1.0;
                        const scalar a3 = gamma_[faceI] - 2.0;
                        const scalar v = (n & U);
                        const scalar vt = meshPhi[faceI] / magSf[faceI];
                        const scalar Lambda = (lc + lv)/magSf[faceI];
    
                        // Calculate A_c * deltaWStar
                        L_rho -= 0.5*magSf[faceI]*( (Lambda - vt)*deltaWStarRho_[own] + (n & deltaWStarRhoU_[own]) );
    
                        L_rhoU.component(0) -= 0.5*magSf[faceI]*
                        (
                            (nx*phi - Ux*v)*deltaWStarRho_[own]
                          + (v - vt - a3*nx*Ux + Lambda)*deltaWStarRhoU_[own].component(0)
                          + (ny*Ux - a2*nx*Uy)*deltaWStarRhoU_[own].component(1)
                          + (nz*Ux - a2*nx*Uz)*deltaWStarRhoU_[own].component(2)
                          + (a2*nx)*deltaWStarRhoE_[own]
                        );

                        L_rhoU.component(1) -= 0.5*magSf[faceI]*
                        (
                            (ny*phi - Uy*v)*deltaWStarRho_[own]
                          + (nx*Uy - a2*ny*Ux)*deltaWStarRhoU_[own].component(0)
                          + (v - vt - a3*ny*Uy + Lambda)*deltaWStarRhoU_[own].component(1)
                          + (nz*Uy - a2*ny*Uz)*deltaWStarRhoU_[own].component(2)
                          + (a2*ny)*deltaWStarRhoE_[own]
                        );
    
                        L_rhoU.component(2) -= 0.5*magSf[faceI]*
                        (
                            (nz*phi - Uz*v)*deltaWStarRho_[own]
                          + (nx*Uz - a2*nz*Ux)*deltaWStarRhoU_[own].component(0)
                          + (ny*Uz - a2*nz*Uy)*deltaWStarRhoU_[own].component(1)
                          + (v - vt - a3*nz*Uz + Lambda)*deltaWStarRhoU_[own].component(2)
                          + (a2*nz)*deltaWStarRhoE_[own]
                        );
    
                        L_rhoE -= 0.5*magSf[faceI]*
                        (
                            v*(phi - a1)*deltaWStarRho_[own]
                          + (nx*a1 - a2*Ux*v)*deltaWStarRhoU_[own].component(0)
                          + (ny*a1 - a2*Uy*v)*deltaWStarRhoU_[own].component(1)
                          + (nz*a1 - a2*Uz*v)*deltaWStarRhoU_[own].component(2)
                          + (gamma_[faceI]*v - vt + Lambda)*deltaWStarRhoE_[own]
                        );
                    }
                }
            }
    
            // Perform calculation for deltaWStar
            deltaWStarRho_[cellI]  = ( -R_rho[cellI]*V[cellI] -  L_rho) / D_[cellI];
            deltaWStarRhoU_[cellI] = (-R_rhoU[cellI]*V[cellI] - L_rhoU) / D_[cellI];
            deltaWStarRhoE_[cellI] = (-R_rhoE[cellI]*V[cellI] - L_rhoE) / D_[cellI];
        }
    
        //
        // Backward sweep
        //

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
    
                        // Arithmetic average
                        const vector Uf = 0.5*(U_[own] + U_[nei]);
                        const scalar cf = 0.5*(c[own] + c[nei]);
                        const scalar nuMaxf = 0.5*(nuMax[own] + nuMax[nei]);
    
                        // Calculate spectral radii
                        const scalar lc = mag( (Uf & Sf[faceI]) - meshPhi[faceI] )
                          + cf*magSf[faceI];
                        const scalar lv = sqr(magSf[faceI]) / dvol * nuMaxf;
    
                        // Face values
                        const scalar rho = rho_[nei] + deltaWRho_[nei];
                        const vector rhoU = rhoU_[nei] + deltaWRhoU_[nei];
                        const scalar rhoE = rhoE_[nei] + deltaWRhoE_[nei];
    
                        const vector U = rhoU / rho;
                        const vector n = Sf[faceI] / magSf[faceI];

                        const scalar Ux = U.component(0);
                        const scalar Uy = U.component(1);
                        const scalar Uz = U.component(2);
                        const scalar nx = n.component(0);
                        const scalar ny = n.component(1);
                        const scalar nz = n.component(2);
    
                        // coefficients
                        const scalar phi = 0.5*(gamma_[faceI] - 1.0)*magSqr(U);
                        const scalar a1 = gamma_[faceI] * rhoE / rho - phi;
                        const scalar a2 = gamma_[faceI] - 1.0;
                        const scalar a3 = gamma_[faceI] - 2.0;
                        const scalar vt = meshPhi[faceI] / magSf[faceI];
                        const scalar v = (n & U);
    
                        // Negative lambda due to upper diagonal
                        const scalar Lambda = -(lc + lv)/magSf[faceI];
    
                        // Calculate A_c * deltaW
                        U_rho += 0.5*( (Lambda - vt)*deltaWRho_[nei] + (n & deltaWRhoU_[nei]) )*magSf[faceI];
    
                        U_rhoU.component(0) += 0.5*magSf[faceI]*
                        (
                            (nx*phi - Ux*v)*deltaWRho_[nei]
                          + (v - vt - a3*nx*Ux + Lambda)*deltaWRhoU_[nei].component(0)
                          + (ny*Ux - a2*nx*Uy)*deltaWRhoU_[nei].component(1)
                          + (nz*Ux - a2*nx*Uz)*deltaWRhoU_[nei].component(2)
                          + (a2*nx)*deltaWRhoE_[nei]
                        );
    
                        U_rhoU.component(1) += 0.5*magSf[faceI]*
                        (
                            (ny*phi - Uy*v)*deltaWRho_[nei]
                          + (nx*Uy - a2*ny*Ux)*deltaWRhoU_[nei].component(0)
                          + (v - vt - a3*ny*Uy + Lambda)*deltaWRhoU_[nei].component(1)
                          + (nz*Uy - a2*ny*Uz)*deltaWRhoU_[nei].component(2)
                          + (a2*ny)*deltaWRhoE_[nei]
                        );
    
                        U_rhoU.component(2) += 0.5*magSf[faceI]*
                        (
                            (nz*phi - Uz*v)*deltaWRho_[nei]
                          + (nx*Uz - a2*nz*Ux)*deltaWRhoU_[nei].component(0)
                          + (ny*Uz - a2*nz*Uy)*deltaWRhoU_[nei].component(1)
                          + (v - vt - a3*nz*Uz + Lambda)*deltaWRhoU_[nei].component(2)
                          + (a2*nz)*deltaWRhoE_[nei]
                        );
    
                        U_rhoE += 0.5*magSf[faceI]*
                        (
                            v*(phi - a1)*deltaWRho_[nei]
                          + (nx*a1 - a2*Ux*v)*deltaWRhoU_[nei].component(0)
                          + (ny*a1 - a2*Uy*v)*deltaWRhoU_[nei].component(1)
                          + (nz*a1 - a2*Uz*v)*deltaWRhoU_[nei].component(2)
                          + (gamma_[faceI]*v - vt + Lambda)*deltaWRhoE_[nei]
                        );
                    }
                }
            }
    
            // Backwards sweep
            deltaWRho_[cellI] = deltaWStarRho_[cellI] - U_rho/D_[cellI];
            deltaWRhoU_[cellI] = deltaWStarRhoU_[cellI] - U_rhoU/D_[cellI];
            deltaWRhoE_[cellI] = deltaWStarRhoE_[cellI] - U_rhoE/D_[cellI];
        }
    
        // Update conserative fields
        rho_ += deltaWRho_;
        rhoU_ += deltaWRhoU_;
        rhoE_ += deltaWRhoE_;

        // Final residuals only on last iteration
        if ( (iter + 1) == nIter_)
        {
            scalar finalRes_rho  = fvc::domainIntegrate( mag(deltaWRho_) / dt ).value();
            scalar finalRes_rhoU = fvc::domainIntegrate( mag(deltaWRhoU_) / dt ).value();
            scalar finalRes_rhoE = fvc::domainIntegrate( mag(deltaWRhoE_) / dt ).value();
        
            Info<< "LUSGS:  Solving for rho,  " 
                << "Initial residual = 1, "
                << "Final residual = " << finalRes_rho/initialRes_rho << ", No Iterations " << nIter_ << nl;
            Info<< "LUSGS:  Solving for rhoU, " 
                << "Initial residual = 1, "
                << "Final residual = " << finalRes_rhoU/initialRes_rhoU << ", No Iterations " << nIter_ << nl;
            Info<< "LUSGS:  Solving for rhoE, "
                << "Initial residual = 1, "
                << "Final residual = " << finalRes_rhoE/initialRes_rhoE << ", No Iterations " << nIter_ << nl;
        }

        // Update primitive fields
        updatePrimitiveFields();
    }
}


// ************************************************************************* //
