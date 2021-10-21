/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) Dr. Martin Heinrich
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
    Foam::functionObjects::performance

Description
    Calculates the performance of a turbomachine between two patches. Resulting
    data consists of (corrected) mass flow rate at first patch, total pressure
    ratio and isentropic efficiency.

    Member function performance::write() calculates the performance and writes
    them into the file \<timeDir\>/performance.dat.

    Example of function object specification:
    \verbatim
    performance1
    {
        type        performance;
        libs        ("libperformance.so");
        
        log         yes;
        patch1      inlet;
        patch2      outlet;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: performance  | yes         |
        log          | Write performance data to standard output | no | no
        patch1       | Patch name 1            | yes |
        patch2       | Patch name 2            | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        corrected    | Optionally correct inlet mass flow rate | no | no
        pRef         | Reference pressure for correction | no | 101325
        scale        | Optionally scale mass flow rate | no | 1
    \endtable

Note
  - If the \c corrected bool is set to true, the mass flow rate at patch 1
    is corrected with respect to the (optional) reference pressure pRef.
    This might be helpful if the inlet pressure is not standard condition.
  - The mass flow rate at patch 1 can be scaled based on the scaling factor
    \c scale. This is useful if only a blade passage of a compressor wheel
    is simulated.

SourceFiles
    performance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_performance_H
#define functionObjects_performance_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class performance Declaration
\*---------------------------------------------------------------------------*/

class performance
:
    public fvMeshFunctionObject,
    public logFiles
{

protected:

    // Protected data

        // Read from dictionary

            //- Inlet and outlet patch name
            label patchID1_;
            label patchID2_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Calculate corrected mass flow rate and reference pressure
            bool corrected_;
            scalar pRef_;

            //- Scale factor used for periodic domains
            label scale_;

            //- Isentropic expansion factor
            scalar gamma_;

            //- Mass flow rate, total pressure ratio and isentropic efficiency
            scalar mDot_;
            scalar Pi_;
            scalar eta_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Calculate the performance 
        void calcPerformance();

public:

    //- Runtime type information
    TypeName("performance");


    // Constructors

        //- Construct from Time and dictionary
        performance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        performance(const performance&) = delete;


    //- Destructor
    virtual ~performance();


    // Member Functions

        //- Read the performance data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the performance
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const performance&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
