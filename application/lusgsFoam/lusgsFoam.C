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
#include "psiThermo.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"

#include "LUSGS.H"
#include "inviscidFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver using the lower-upper"
        " symmetric Gauss Seidel implicit time integration schemes"
        " and Kurvanov / Tadmor inviscid flux schemes."
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        #include "setDeltaT.H"
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.update();
        
        // Update fluxes
        flux.update();

        // Get acoustic Courant number
        CoNum = flux.CoNum();

        // Perform forward & backward sweep of LU-SGS scheme
        lusgs.sweep();

        // Update the mass flux            
        phi = flux.phi();

        // Solve governing equations of turbulence model
        turbulence->correct();

        // Write time step
        runTime.write();

      	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "\nEnd\n" << endl;
}

// ************************************************************************* //
