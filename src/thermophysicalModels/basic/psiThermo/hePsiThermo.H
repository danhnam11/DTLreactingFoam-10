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
    Foam::hePsiThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    hePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiThermo_H
#define hePsiThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class hePsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables, for sutherland only
        void calculate();
        // Nam 
        //- Initialize the thermo variables with full DTM/FTM
        void initialize();
        //- Calculate the thermo variables with coTHERM
        void calculateUsingCoTHERM();
        //- Calculate the thermo variables with coTHERM only T
        void calculateUsingCoTHERMOnlyT();
        //- Calculate the thermo variables, only for preprocessing FTM
        void calculateTransportPreProcessing();
        // Nam

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


public:

    //- Runtime type information
    TypeName("hePsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        hePsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        hePsiThermo(const hePsiThermo<BasicPsiThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~hePsiThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


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




    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const hePsiThermo<BasicPsiThermo, MixtureType>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hePsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
