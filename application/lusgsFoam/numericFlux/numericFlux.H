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

SourceFiles
    numericFlux

Description
    Central-upwind schemes of Kurganov and Tadmor for flux calculation. Uses
    primitive variables as input and gives back conservative numerical fluxes.
 
    A. Kurganov, E. Tadmor. "New High-Resolution Central Schemes for Nonlinear
    Conservation Laws and Convection-Diffusion Equations."
    Journal of Computational Physics 160 (2000): 241-282.

    Kurganov, S. Noelle, G. Petrova. "Semi-discrete central-upwind schemes for
    hyperbolic conservation laws and Hamilton-Jacobi equations."
    SIAM Journal on Scientific Computing 23 (2001): 707-740. 

\*---------------------------------------------------------------------------*/

#ifndef numericFlux_H
#define numericFlux_H

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class numericFlux Declaration
\*---------------------------------------------------------------------------*/

class numericFlux
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        // Reference to primitive fields

            //- Pressure
            const volScalarField& p_;

            //- Velocity
            const volVectorField& U_;

            //- Density
            const volScalarField& rho_;

        //- Reference to the thermophysical model
        const psiThermo& thermo_;

        //- Reference to turbulence model
        const compressible::turbulenceModel& turbulence_;

        // Fluxes
        
            surfaceScalarField phi_;
            surfaceVectorField phiUp_;
            surfaceScalarField phiEp_;
            surfaceScalarField amaxSf_;

            //- Positive & negative direction
            const surfaceScalarField pos_;
            const surfaceScalarField neg_;

        // Reynolds stress tensor
        volTensorField tauMC_;

        // Dissipation
        surfaceScalarField sigmaDotU_;
                        
        // Constant
        const dimensionedScalar vZero_;
        
        // Flux scheme
        word fluxScheme_;
    


public:

    // Constructors

        //- Construct from components
        numericFlux
        (
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& rho,
            const psiThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );

    // Member Functions

        // Return variables
 
            // Mass flux
            const surfaceScalarField& phi()
            {
                return phi_;
            }

            // Momentum flux
            const surfaceVectorField& phiUp() const
            {
                return phiUp_;
            }

            // Energy flux
            const surfaceScalarField& phiEp() const
            {
                return phiEp_;
            }

            // Viscous stress tensor
            const volTensorField& tauMC() const
            {
                return tauMC_;
            }

            // Dissipation
            const surfaceScalarField& sigmaDotU() const
            {
                return sigmaDotU_;
            }

        //- Update fluxes
        void update();

        //- Print Courant number
        void courantNo();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
