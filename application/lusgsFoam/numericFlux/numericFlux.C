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

#include "numericFlux.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::numericFlux::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& rho,
    const psiThermo& thermophysicalModel,
    const compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    rho_(rho),
    thermo_(thermophysicalModel),
    turbulence_(turbulenceModel),
    phi_
    (
        IOobject
        (
            "flux::phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("flux::phi", dimMass/dimTime, scalar(0.0))
    ),
    amaxSf_
    (
        IOobject
        (
            "flux::amaxSf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("flux::amaxSf", dimVolume/dimTime, scalar(0.0))
    ),
    pos_
    (
        IOobject
        (
            "flux::pos",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("pos", dimless, 1.0)
    ),
    neg_
    (
        IOobject
        (
            "flux::neg",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("neg", dimless, -1.0)
    ),
    resRho_
    (
        IOobject
        (
            "flux::resRho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("flux::resRho", dimMass/dimTime/dimVolume, scalar(0.0))
    ),
    resRhoU_
    (
        IOobject
        (
            "flux::resRhoU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("flux::resRhoU", dimMass/dimTime/dimTime/dimArea, vector(0,0,0))
    ),
    resRhoE_
    (
        IOobject
        (
            "flux::resRhoE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("flux::resRhoE", dimMass/dimLength/dimTime/dimTime/dimTime, scalar(0.0))
    ),
    vZero_
    (
        dimensionedScalar("vZero", dimVolume/dimTime, 0.0)
    )
{
    fluxScheme_ = "Kurganov";
    if (mesh_.schemesDict().readIfPresent("fluxScheme", fluxScheme_))
    {
        if ((fluxScheme_ == "Tadmor") || (fluxScheme_ == "Kurganov"))
        {
            Info<< "fluxScheme: " << fluxScheme_ << endl;
        }
        else
        {
            FatalErrorInFunction
                << "fluxScheme: " << fluxScheme_
                << " is not a valid choice. "
                << "Options are: Tadmor, Kurganov"
                << abort(FatalError);
        }
    }

    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::numericFlux::update()
{
    // Reconstruct velocity
    surfaceVectorField U_pos(interpolate(U_, pos_)); 
    surfaceVectorField U_neg(interpolate(U_, neg_));
    
    // Reconstruct pressure
    surfaceScalarField p_pos(interpolate(p_, pos_));
    surfaceScalarField p_neg(interpolate(p_, neg_));
    
    // Reconstruct density
    surfaceScalarField rho_pos(interpolate(rho_, pos_));
    surfaceScalarField rho_neg(interpolate(rho_, neg_));
    
    // Reconstruct internal energy
    surfaceScalarField e_pos(interpolate(thermo_.he(), pos_, p_.name()));
    surfaceScalarField e_neg(interpolate(thermo_.he(), neg_, p_.name()));
    
    // Volumetric flux
    surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh_.Sf());
    surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh_.Sf());

    // Make fluxes relative to mesh-motion
    if (mesh_.moving())
    {
        phiv_pos -= mesh_.phi();
        phiv_neg -= mesh_.phi();
    }
   
    // Note: extracted out the orientation so becomes unoriented
    phiv_pos.setOriented(false);
    phiv_neg.setOriented(false);

    volScalarField c = sqrt(thermo_.Cp()/thermo_.Cv() / thermo_.psi());

    // Sonic flux
    surfaceScalarField cSf_pos
    (
        "cSf_pos",
        interpolate(c, pos_, rho_.name())*mesh_.magSf()
    );
    surfaceScalarField cSf_neg
    (
        "cSf_neg",
        interpolate(c, neg_, rho_.name())*mesh_.magSf()
    );

    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), vZero_)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), vZero_)
    );

    surfaceScalarField a_pos("a_pos", ap/(ap - am));
    
    amaxSf_ = max(mag(am), mag(ap));

    surfaceScalarField aSf("aSf", am*a_pos);

    if (fluxScheme_ == "Tadmor")
    {
        aSf = -0.5*amaxSf_;
        a_pos = 0.5;
    }
    surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

    phiv_pos *= a_pos;
    phiv_neg *= a_neg;
    
    surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
    surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);
    
    // Reuse amaxSf for the maximum positive and negative fluxes
    // estimated by the central scheme
    amaxSf_ = max(mag(aphiv_pos), mag(aphiv_neg));
    
    surfaceVectorField phiU(aphiv_pos*rho_pos*U_pos + aphiv_neg*rho_neg*U_neg);
    // Note: reassembled orientation from the pos and neg parts so becomes
    // oriented
    phiU.setOriented(true);
    
    // Update fluxes
    phi_ = aphiv_pos*rho_pos + aphiv_neg*rho_neg;
    phiU += (a_pos*p_pos + a_neg*p_neg)*mesh_.Sf();
    surfaceScalarField phiE = aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg;

    phi_.setOriented(true);

    // Make flux for pressure-work absolute
    if (mesh_.moving())
    {
        surfaceScalarField phia(a_pos*p_pos + a_neg*p_neg);
        phia.setOriented(true);

        phiE += mesh_.phi()*phia;
    }

    // Viscous and Reynolds stress tensor
    const volSymmTensorField tau
    (
        "tau",
      - turbulence_.devRhoReff()
      - ((2.0/3.0)*I)*rho_*turbulence_.k()
    );

    // Assemble residuals for governing equations
    resRho_ =
        fvc::div(phi_);
    resRhoU_ =
        fvc::div(phiU)
      - fvc::div(tau);
    resRhoE_ = 
        fvc::div(phiE)
      - fvc::div(tau & U_, "div(tau&U)")
      - fvc::laplacian
        (
            turbulence_.mu() + 0.6*turbulence_.mut(),
            turbulence_.k()
        )
      - fvc::laplacian(turbulence_.alphaEff(), thermo_.he());
}


void Foam::numericFlux::courantNo()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    if (mesh_.nInternalFaces())
    {
        scalarField sumAmaxSf(fvc::surfaceSum(amaxSf_)().primitiveField());

        CoNum = 0.5*gMax(sumAmaxSf/mesh_.V().field())*mesh_.time().deltaTValue();

        meanCoNum =
            0.5*(gSum(sumAmaxSf)/gSum(mesh_.V().field()))*mesh_.time().deltaTValue();
    }

    Info<< "Acoustic Courant mean:" << meanCoNum
        << " max: " << CoNum << endl;
}


// ************************************************************************* //
