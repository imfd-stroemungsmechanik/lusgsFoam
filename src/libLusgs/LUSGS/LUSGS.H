/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 Dr. Martin Heinrich
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
    Foam::LUSGS

Description
    Implicit time advancement algorithm based on the Lower-Upper Symmetric-
    Gauss-Seidel algorithm.

\*---------------------------------------------------------------------------*/

#ifndef LUSGS_H
#define LUSGS_H

#include "fvCFD.H"
#include "psiThermo.H"
#include "numericFlux.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"

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
        
        //- Reference to thermophysical transport model
        const fluidThermophysicalTransportModel& thermophysicalTransport_;

        //- Reference to turbulence model
        const compressible::momentumTransportModel& turbulence_;

        //- Reference to MRFZones
        const MRFZoneList& MRF_;

        //- LUSGS solver settings
        
            //- Over-relaxation factor
            const scalar omega_;

            //- Switch for inviscid flows
            bool inviscid_;

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
            const fluidThermophysicalTransportModel& thermophysicalTransport,
            const compressible::momentumTransportModel& turbulence,
            const IOMRFZoneList& MRF
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
