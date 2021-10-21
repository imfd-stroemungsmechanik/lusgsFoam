/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM:
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of OpenFOAM nor foam-extend.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
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

SourceFiles
    numericFlux.H
    numericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef numericFlux_H
#define numericFlux_H

#include "fvCFD.H"
#include "psiThermo.H"
#include "MRFZoneList.H"
#include "directionInterpolate.H"
#include "dynamicMomentumTransportModel.H"

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
    
            //- Static pressure
            const volScalarField& p_;
    
            //- Velocity
            const volVectorField& U_;
    
            //- Density
            const volScalarField& rho_;
    
        //- Reference to the thermophysicalModel
        const psiThermo& thermo_;

        //- Reference to turbulence model
        const compressible::momentumTransportModel& turbulence_;
    
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
            const compressible::momentumTransportModel& turbulence,
            const MRFZoneList& MRF
        );
    
    
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
    
        //- Update fluxes 
        void update();

        //- Return Courant number
        scalar CoNum()
        {
            return CoNum_;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
