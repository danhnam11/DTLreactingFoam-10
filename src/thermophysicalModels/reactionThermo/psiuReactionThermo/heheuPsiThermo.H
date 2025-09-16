/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::heheuReactionThermo

Description
    Foam::heheuReactionThermo

SourceFiles
    heheuReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heheuPsiThermo_H
#define heheuPsiThermo_H

#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class heheuPsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heheuPsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Data

        //- Unburnt temperature
        volScalarField Tu_;

        //- Unburnt energy
        volScalarField heu_;


        // Nam - for diffusivity models in DTM 
        //- Store list of diffusion coefficients
        PtrList<volScalarField> Dimix_;

        //- Store list of thermal diffusion coefficients
        PtrList<volScalarField> DimixT_;

        //- Store list of binary diffusion coefficients
        List<PtrList<volScalarField>> Dij_;

        //- Store List of energy field of indivisual species [J/kg] 
        PtrList<volScalarField> heList_;

        //- Store molecular weight of mixture
        volScalarField Wmix_;

        //- Store List of raw thermophysical data of all individual species 
        const PtrList<typename MixtureType::thermoType>& specieThermos_;
        // Nam - end

    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


public:

    //- Runtime type information
    TypeName("heheuPsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        heheuPsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        heheuPsiThermo
        (
            const heheuPsiThermo<BasicPsiThermo, MixtureType>&
        ) = delete;


    //- Destructor
    virtual ~heheuPsiThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables.

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu()
            {
                return heu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const
            {
                return heu_;
            }

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const
            {
                return Tu_;
            }


        // Fields derived from thermodynamic state variables

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const label patchi
            ) const;


            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const;


        // Access to transport variables

            //- Dynamic viscosity of unburnt gas [kg/m/s]
            virtual tmp<volScalarField> muu() const;

            //- Dynamic viscosity of burnt gas [kg/m/s]
            virtual tmp<volScalarField> mub() const;


        // Nam - for diffusivity models in DTM 
        // - Diffusion coefficient of specie ith in the mixture [m^2/s] - Nam
        virtual tmp<volScalarField> Dimix( const label speciei) const
        {
            return Dimix_[speciei];
        }

        // - Diffusion coefficient of specie ith in the mixture for patch [m^2/s] - Nam
        virtual tmp<scalarField> Dimix( const label speciei, const label patchi) const
        {
            return Dimix_[speciei].boundaryField()[patchi];
        }

        // - Thermal diffusion coefficient of specie ith in the mixture [m^2/s] - Nam
        virtual tmp<volScalarField> DimixT( const label speciei) const
        {
            return DimixT_[speciei];
        }

        // - Thermal diffusion coefficient of specie ith in the mixture for patch [m^2/s] - Nam
        virtual tmp<scalarField> DimixT( const label speciei, const label patchi) const
        {
            return DimixT_[speciei].boundaryField()[patchi];
        }

        // - Binary diffusion coefficient of species ith and species jth in the mixture [m^2/s] - JH
        virtual tmp<volScalarField> Dij( const label speciei, const label speciej) const
        {
            return Dij_[speciei][speciej];
        }

        // - Binary diffusion coefficient of species ith and species jth in the mixture for patch [m^2/s] - JH
        virtual tmp<scalarField> Dij( const label speciei, const label speciej,  const label patchi) const
        {
            return Dij_[speciei][speciej].boundaryField()[patchi];
        }

        //- Return he[J/kg] of indivisual speciei
        virtual tmp<volScalarField> hei(label speciei) const
        { 
            return heList_[speciei];
        }  

        //- Return he[J/kg] of indivisual speciei
        virtual tmp<scalarField> hei(label speciei, const label patchi) const
        {
            return heList_[speciei].boundaryField()[patchi];
        }

        // return molecular weight of mixture 
        virtual tmp<volScalarField> Wmix() const
        {
            return Wmix_; 
        }
        // Nam 


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "heheuPsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
