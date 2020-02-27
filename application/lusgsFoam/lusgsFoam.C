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
#include "turbulentFluidThermoModel.H"

#include "numericFlux.H"
#include "LUSGS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"
        flux.courantNo();

        // Iteration counter
        int nIter(0);

        // Start inner correctors
        for (int iter=0; iter < maxIter; iter++)
        {
            nIter++;

            // Update fluxes
            flux.update();

            // Viscous and Reynolds stress tensor
            const volSymmTensorField tau
            (
                "tau",
              - turbulence->devRhoReff()
              - ((2.0/3.0)*I)*rho*turbulence->k()
            );

            // Assemble residuals
            volScalarField R_rho
            (
                fvc::div(flux.phi())
            );
            volVectorField R_rhoU
            (
                fvc::div(flux.phiU())
              - fvc::div(tau)
            );
            volScalarField R_rhoE
            (
                fvc::div(flux.phiE())
              - fvc::div(tau & U, "div(tau&U)")
              - fvc::laplacian
                (
                    (turbulence->mu()+0.6*turbulence->mut()),
                    turbulence->k()
                )
              - fvc::laplacian(turbulence->alphaEff(), e)
            );

            // Perform forward & backward sweep of LU-SGS scheme
            lusgs.update(R_rho, R_rhoU, R_rhoE);

            // Store initial and current residual
            if (iter == 0)
            {
	            initRho = lusgs.resRho();
	            initRhoU = lusgs.resRhoU();
	            initRhoE = lusgs.resRhoE();
	        }

            // Current residual
            curRho = lusgs.resRho();
            curRhoU = lusgs.resRhoU();
            curRhoE = lusgs.resRhoE();

            // Break correctors if relative residual is reached
            if
            (
                (curRho/initRho) < relTol
             && (curRhoU/initRhoU) < relTol
             && (curRhoE/initRhoE) < relTol
            )
            {
                break;
            }
        }

        Info<< "LUSGS:  Solving for rho,  " 
            << "Initial residual = " << initRho << ", "
            << "Final residual = " << curRho << ", No Iterations " << nIter << nl;
        Info<< "LUSGS:  Solving for rhoU,  " 
            << "Initial residual = " << initRhoU << ", "
            << "Final residual = " << curRhoU << ", No Iterations " << nIter << nl;
        Info<< "LUSGS:  Solving for rhoE,  " 
            << "Initial residual = " << initRhoE << ", "
            << "Final residual = " << curRhoE << ", No Iterations " << nIter << nl;

        // Update mass flux
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
