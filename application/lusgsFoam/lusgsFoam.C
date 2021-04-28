/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    lusgsFoam

Description
    Density-based solver for high-speed, transient flows. Uses implicit
    time advancement for large Courant-numbers based on the Lower-Upper
    Symmetric Gauss-Seidel (LU-SGS) scheme.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "motionSolver.H"

#include "numericFlux.H"
#include "LUSGS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.update();

        #include "compressibleCourantNo.H"
        flux.courantNo();

        Info << "Flux calculation... ";
        flux.update();

        // Make flux accessible for other parts of openfoam
        phi = flux.phi();

        // Perform forward & backward sweep of LU-SGS scheme
        lusgs.sweep();

        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "\nEnd\n" << endl;
}

// ************************************************************************* //
