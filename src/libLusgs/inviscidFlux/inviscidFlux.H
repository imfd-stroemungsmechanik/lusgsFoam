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
    inviscidFlux

Description
    Inviscid flux class, computes inviscid numerical fluxes based on the
    central-upwind schemes of Kurganov and Tadmor. Uses primitive variables
    as input and gives back conservative numerical fluxes.

    Kurganov A, Tadmor E.: New high-resolution central schemes for nonlinear
    conservation laws and convection–diffusion equations.
    Journal of Computational Physics 2001; 160:241–282.

\*---------------------------------------------------------------------------*/

#ifndef inviscidFlux_H
#define inviscidFlux_H

#include "fvCFD.H"
#include "psiThermo.H"
#include "MRFZoneList.H"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inviscidFlux Declaration
\*---------------------------------------------------------------------------*/

class inviscidFlux
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        // Reference to primitive fields

            //- static Pressure
            const volScalarField& p_;

            //- Velocity
            const  volVectorField& U_;

            //- static density
            const  volScalarField& rho_;

        //- Reference to the thermophysicalModel
        const psiThermo& thermo_;

        //- Reference to turbulence model
        const compressible::turbulenceModel& turbulence_;
    
        //- Reference to MRFZones
        const MRFZoneList& MRF_;

        //- Isentropic expansion factor
        const surfaceScalarField gamma_;

        // Fluxes
    
            //- Mass flux
            surfaceScalarField phi_;
    
            //- Momentum flux
            surfaceVectorField phiUp_;
    
            //- Energy flux
            surfaceScalarField phiEp_;

            //- Positive & negative direction
            const surfaceScalarField pos_;
            const surfaceScalarField neg_;

        // Reynolds stress tensor
        volTensorField tauMC_;

        // Dissipation
        surfaceScalarField sigmaDotU_;

        // Constant
        const dimensionedScalar vZero_;

        //- Maximum and mean acoustic Courant number
        scalar CoNum_;
        scalar meanCoNum_;
        
        // Flux scheme
        const word fluxScheme_;

public:

    // Constructors

        //- Construct from components
        inviscidFlux
        (
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& rho,
            const psiThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const MRFZoneList& MRF
        );


    // Destructor - default

    // Member Functions
    
        // Return variables
    
            //- Return mass flux
            const surfaceScalarField& phi() const
            {
                return phi_;
            }
    
            //- Return momentum flux
            const surfaceVectorField& phiUp() const
            {
                return phiUp_;
            }
    
            //- Return energy flux
            const surfaceScalarField& phiEp() const
            {
                return phiEp_;
            }
    
            //- Return turbulent shear stress
            const volTensorField& tauMC() const
            {
                return tauMC_;
            }
    
            //- Return dissipation
            const surfaceScalarField& sigmaDotU() const
            {
                return sigmaDotU_;
            }
    
        //- Update convective fluxes 
        void update();

        //- Return Courant number
        scalar CoNum();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
