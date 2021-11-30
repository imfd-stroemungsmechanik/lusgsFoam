/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::lossCoefficient

Description
    Calculates the lossCoefficient of a turbomachine between two patches. Resulting
    data consists of (corrected) mass flow rate at first patch, total pressure
    ratio and isentropic efficiency.

    Member function lossCoefficient::write() calculates the lossCoefficient and writes
    them into the file \<timeDir\>/lossCoefficient.dat.

    Example of function object specification:
    \verbatim
    lossCoefficient1
    {
        type        lossCoefficient;
        libs        ("liblossCoefficient.so");
        
        log         yes;
        patch1      inlet;
        patch2      outlet;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: lossCoefficient  | yes         |
        log          | Write lossCoefficient data to standard output | no | no
        patch1       | Patch name 1            | yes |
        patch2       | Patch name 2            | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
    \endtable

SourceFiles
    lossCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef lossCoefficient_H
#define lossCoefficient_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "fluidThermo.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class lossCoefficient Declaration
\*---------------------------------------------------------------------------*/

class lossCoefficient
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        // File streams
        autoPtr<OFstream> lossCoefficientFilePtr_;

        //- Inlet and outlet patch name
        label patchID1_;
        label patchID2_;

        //- Name of pressure field
        word pName_;

        //- Name of velocity field
        word UName_;

        // Loss coefficient
        scalar lossCoefficient_;

        //- Isentropic expansion factor
        scalar gamma_;

    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Calculate the lossCoefficient 
        void calclossCoefficient();

        //- Disallow default bitwise copy construct
        lossCoefficient(const lossCoefficient&);

        //- Disallow default bitwise assignment
        void operator=(const lossCoefficient&);


public:

    //- Runtime type information
    TypeName("lossCoefficient");


    // Constructors

        //- Construct from Time and dictionary
        lossCoefficient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
        
        //- Construct from objectRegistry and dictionary
        lossCoefficient
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    //- Destructor
    virtual ~lossCoefficient();


    // Member Functions

        //- Read the lossCoefficient data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the acoustic data
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
