/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "performance.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(performance, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        performance,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::performance::writeFileHeader(const label i)
{
    writeHeader(file(), "performance");

    writeCommented(file(), "Time");
    writeTabbed(file(), "(corrected) mass flow rate");
    writeTabbed(file(), "total pressure ratio");
    writeTabbed(file(), "isentropic efficiency");
    file() << endl;
}


void Foam::functionObjects::performance::calcPerformance()
{
    mDot_ = 0.0;
    Pi_ = 0.0;
    eta_ = 0.0;

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    const volScalarField& T = obr_.lookupObject<volScalarField>("T");
    const surfaceScalarField& phi = obr_.lookupObject<surfaceScalarField>("phi");

    const volScalarField& psi = obr_.lookupObject<volScalarField>("thermo:psi");
    const volScalarField Ma = mag(U)/sqrt(gamma_/psi);

    const volScalarField pTot
    (
        p * pow( ( 1.0 + (gamma_ - 1.0) / 2.0 * Ma*Ma), (gamma_ / (gamma_ - 1.0)))
    );

    const volScalarField Ttot
    (
        T * ( 1.0 + 0.5 * (gamma_ - 1.0) * Ma*Ma)
    );

    // Mass flow rate at the inlet and outlet
    mDot_ = mag(gSum(phi.boundaryField()[patchID1_]));
    const scalar mDotOut = mag(gSum(phi.boundaryField()[patchID2_]));

    // mass flow averaged total pressure
    scalar p01 = mag(gSum( pTot.boundaryField()[patchID1_] * phi.boundaryField()[patchID1_] ) ) / max(mDot_, VSMALL);
    scalar p02 = mag(gSum( pTot.boundaryField()[patchID2_] * phi.boundaryField()[patchID2_] ) ) / max(mDotOut, VSMALL);

    // mass flow averaged total pressure
    scalar T01 = mag(gSum( Ttot.boundaryField()[patchID1_] * phi.boundaryField()[patchID1_] ) ) / max(mDot_, VSMALL);
    scalar T02 = mag(gSum( Ttot.boundaryField()[patchID2_] * phi.boundaryField()[patchID2_] ) ) / max(mDotOut, VSMALL);

    if (corrected_)
    {
        mDot_ *= pRef_ / p01;
    }

    // Account for scale
    mDot_ *= scale_;

    Pi_ = p02 / p01;
    eta_ = max
    (
        0.0,
        min
        (
            (Foam::pow( p02 / p01, (gamma_ - 1.0) / gamma_ ) - 1.0) / ( T02 / T01 - 1.0),
            1.0
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::performance::performance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    patchID1_(Zero),
    patchID2_(Zero),
    pName_(word::null),
    UName_(word::null),
    corrected_(false),
    pRef_(Zero),
    scale_(Zero),
    gamma_(Zero),
    mDot_(Zero),
    Pi_(Zero),
    eta_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::performance::~performance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::performance::read(const dictionary& dict)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    fvMeshFunctionObject::read(dict);

    resetName(typeName);

    // Inlet and outlet patch names
    word patchName1 = word(dict.lookup("patch1"));
    word patchName2 = word(dict.lookup("patch2"));

    patchID1_ = mesh.boundaryMesh().findPatchID(patchName1);
    patchID2_ = mesh.boundaryMesh().findPatchID(patchName2);

    // Optional U, p entries
    pName_ = dict.lookupOrDefault<word>("pName", "p");
    UName_ = dict.lookupOrDefault<word>("UName", "U");

    // Corrected mass flow rate calculation
    corrected_ = dict.lookupOrDefault<bool>("corrected", false);
    if (corrected_)
    {
        pRef_ = readScalar(dict.lookup("pRef"));
    }

    // Scale factor
    scale_ = dict.lookupOrDefault<label>("scale", 1);

    // Isentropic expansion factor
    gamma_ = dict.lookupOrDefault<scalar>("gamma", 1.4);

    if
    (
        !obr_.foundObject<volVectorField>(UName_)
     || !obr_.foundObject<volScalarField>(pName_)
    )
    {
        FatalErrorInFunction
            << "Could not find " << UName_ << ", " << pName_
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::performance::execute()
{
    calcPerformance();

    return true;
}


bool Foam::functionObjects::performance::write()
{
    Log << type() << " " << name() << " write:" << nl;

    logFiles::write();

    if (Pstream::master())
    {
        Log << "    (corrected) mass flow rate : " << mDot_ << nl
            << "          total pressure ratio : " << Pi_ << nl
            << "         isentropic efficiency : " << eta_
            << endl;

        writeTime(file());

        file()
            << tab << mDot_
            << tab << tab << Pi_
            << tab << eta_
            << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
