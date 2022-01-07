/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "lossCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(lossCoefficient, 0);

    addToRunTimeSelectionTable(functionObject, lossCoefficient, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::lossCoefficient::createFiles()
{
    if (writeToFile() && !lossCoefficientFilePtr_.valid())
    {
        lossCoefficientFilePtr_ = createFile("lossCoefficientType");
        writeIntegratedHeader("lossCoefficientType", lossCoefficientFilePtr_());
    }
}


void Foam::functionObjects::lossCoefficient::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeCommented(os, "Time");
    writeTabbed(os, "loss coefficient [-]");
    os  << endl;
}


void Foam::functionObjects::lossCoefficient::calclossCoefficient()
{
    lossCoefficient_ = 0.0;

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = obr_.lookupObject<surfaceScalarField>("phi");

    const volScalarField& psi = obr_.lookupObject<volScalarField>("thermo:psi");
    const volScalarField& rho = obr_.lookupObject<volScalarField>("rho");
    const volScalarField Ma = mag(U)/sqrt(gamma_/psi);

    const scalar gM1ByG = (gamma_ - 1.0)/gamma_;

    const volScalarField pTot
    (
        p * pow( ( 1.0 + (gamma_ - 1.0) / 2.0 * Ma*Ma), (gamma_ / (gamma_ - 1.0)))
    );


    // Mass flow rate at the inlet and outlet
    const scalar m1 = mag(gSum(phi.boundaryField()[patchID1_]));
    const scalar m2 = mag(gSum(phi.boundaryField()[patchID2_]));

    if ( (mag(m1) > SMALL) && (mag(m2) > SMALL) )
    {
        // mass flow averaged static pressure
        scalar p1 = mag(gSum( p.boundaryField()[patchID1_] * phi.boundaryField()[patchID1_] ) ) / max(m1, SMALL);
        scalar p2 = mag(gSum( p.boundaryField()[patchID2_] * phi.boundaryField()[patchID2_] ) ) / max(m2, SMALL);

        // mass flow averaged total pressure
        scalar p01 = mag(gSum( pTot.boundaryField()[patchID1_] * phi.boundaryField()[patchID1_] ) ) / max(m1, SMALL);
        scalar p02 = mag(gSum( pTot.boundaryField()[patchID2_] * phi.boundaryField()[patchID2_] ) ) / max(m2, SMALL);

        lossCoefficient_ = 1 - ( 1 - pow(p2 / p02, (gamma_ - 1.0) / gamma_) )
                        / ( 1 - pow(p2 / p01, (gamma_ - 1.0) / gamma_) );

    }
    else
    {
        lossCoefficient_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::lossCoefficient::lossCoefficient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    patchID1_(Zero),
    patchID2_(Zero),
    pName_(word::null),
    UName_(word::null),
    lossCoefficient_(Zero),
    gamma_(Zero)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::functionObjects::lossCoefficient::~lossCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::lossCoefficient::read(const dictionary& dict)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Inlet and outlet patch names
    word patchName1 = word(dict.lookup("patch1"));
    word patchName2 = word(dict.lookup("patch2"));

    patchID1_ = mesh.boundaryMesh().findPatchID(patchName1);
    patchID2_ = mesh.boundaryMesh().findPatchID(patchName2);

    // Optional U, p entries
    pName_ = dict.lookupOrDefault<word>("pName", "p");
    UName_ = dict.lookupOrDefault<word>("UName", "U");

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

    // Isentropic expansion factor
    gamma_ = dict.lookupOrDefault<scalar>("gamma", 1.4);

    return true;
}


bool Foam::functionObjects::lossCoefficient::execute()
{
    calclossCoefficient();

    return true;
}


bool Foam::functionObjects::lossCoefficient::write()
{

    if (Pstream::master())
    {
        createFiles();

        Info<< type() << " output:" << nl
            << "        loss coefficient  : " << lossCoefficient_ << nl
            << endl;

        // File output
        if (writeToFile())
        {
            Ostream& os = lossCoefficientFilePtr_();

            writeCurrentTime(os);
            
            os  << tab << lossCoefficient_ << tab
                << endl;
        }
    }
    return true;
}


// ************************************************************************* //
