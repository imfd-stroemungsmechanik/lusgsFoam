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

Class
    Foam::LUSGS

Description
    Implicit time advancement algorithm based on the Lower-Upper Symmetric-
    Gauss-Seidel algorithm. Exact derivations of invicid flux Jacobians
    are implemented. 

\*---------------------------------------------------------------------------*/

#ifndef LUSGS_H
#define LUSGS_H

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "numericFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LUSGS Declaration
\*---------------------------------------------------------------------------*/

class LUSGS
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to primitive fields

            //- Static Pressure
            volScalarField& p_;

            //- Density
            volScalarField& rho_;

            //- Velocity
            volVectorField& U_;

        //- Reference to numeric flux
        numericFlux& flux_;

        //- Reference to thermophysical model
        psiThermo& thermo_;
        
        //- Reference to turbulence model
        const compressible::turbulenceModel& turbulence_;

        //- LUSGS solver settings
        
            //- Over-relaxation factor
            const scalar omega_;

        //- Isentropic expansion factor
        const surfaceScalarField gamma_;
            
        //- Diagonal matrix
        scalarField D_;
            
        //- Conservative fields
        
            //- Momentum
            volVectorField rhoU_;
            
            //- Energy
            volScalarField rhoE_;

        //- Intermediate fields
        
            //- Forward sweep
            volScalarField deltaWStarRho_;
            volVectorField deltaWStarRhoU_;
            volScalarField deltaWStarRhoE_;

            //- Backward sweep
            volScalarField deltaWRho_;
            volVectorField deltaWRhoU_;
            volScalarField deltaWRhoE_;

    // Private Member Functions
        
        //- Update primitive fields
        void updatePrimitiveFields();

public:

    // Constructors

        //- Construct from components
        LUSGS
        (
            volScalarField& p,
            volScalarField& rho,
            volVectorField& U,
            numericFlux& flux,
            psiThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    // Member Functions

        //- Perform forward & backward sweep
        void sweep();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
