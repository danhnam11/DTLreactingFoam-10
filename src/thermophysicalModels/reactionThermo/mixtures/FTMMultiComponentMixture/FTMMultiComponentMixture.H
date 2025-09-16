/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::FTMMultiComponentMixture

Description

    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to the thermodynamic coefficients.

    Transport properties are based on Standard Kinetic Theory model. 
    The mixture is calculated based on the modification of the Wilke semi-empirical formulas. 

    Mass Diffusivity is based on mixture averaged model in which binary diffusion 
    coefficients are obtained from Standard Kinetic Theory model.

    It has been validated against NIST, CHEMKIN, Cantera, laminarSMOKE
 
    For details of equations in the model, refer to TRANSPORT program 
    in CHEMKIN II @ R. J. Kee et al. CHEMKIN collection, Release 3.6, 
    Reaction Design, Inc., San Diego, CA (2001). 

    Output: 
        Dynamic Viscosity [kg/m.s]
        Thermal Conductivity [W/m.K]    
        Mass Diffusivity [m2/s]    

    developed by:
    Nam Danh Nguyen and Jae Hun Lee
    Advisor: Prof. Chun Sang Yoo 
    Clean combustion & energy research (CCER) Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    FTMMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef FTMMultiComponentMixture_H
#define FTMMultiComponentMixture_H

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class FTMMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class FTMMultiComponentMixture
:
    public multiComponentMixture<ThermoType>
{

public:

    //- Mixing type for thermodynamic properties
    typedef ThermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef ThermoType transportMixtureType;


private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoMixtureType mixture_;


        //- For mixture DTM - Nam
        //- List of species molecular weight
        List<scalar> ListW_;

        List<List<scalar>> muCoeffsMk_;
        List<List<scalar>> kappaCoeffsMk_;
        List<List<List<scalar>>> DijCoeffsMk_;
        // 


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        FTMMultiComponentMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        FTMMultiComponentMixture
        (
            const FTMMultiComponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~FTMMultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "FTMMultiComponentMixture<" + ThermoType::typeName() + '>';
        }

        const thermoMixtureType& cellThermoMixture(const label celli) const;

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label patchi,
            const label facei
        ) const;


        const transportMixtureType& cellTransportMixture
        (
            const label celli
        ) const
        {
            return cellThermoMixture(celli);
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return patchFaceThermoMixture(patchi, facei);
        }

        const transportMixtureType& cellTransportMixture
        (
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FTMMultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
