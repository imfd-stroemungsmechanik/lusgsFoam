Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
thermo.validate(args.executable(), "e");

volScalarField& p = thermo.p();
volScalarField& e = thermo.he();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Create mass flux alias for easier coupling with other code components
surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rho*U));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

#include "createMRF.H"

// Flux calculation
inviscidFlux flux
(
    p,
    U,
    rho,
    thermo,
    turbulence,
    MRF
);

LUSGS lusgs
(
    p,
    U,
    rho,
    flux,
    thermo,
    turbulence,
    MRF
);
