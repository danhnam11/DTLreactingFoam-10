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
    Foam::fluidReactionThermo

Description
    Base-class for multi-component fluid thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    fluidReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidReactionThermo_H
#define fluidReactionThermo_H

#include "fluidThermo.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluidReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidReactionThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("fluidReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~fluidReactionThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;



        // New functions added for diffusivity models - Nam
        // - Diffusion coefficient of specie ith in the mixture [m ^2/ s]
        virtual tmp<volScalarField> Dimix(const label speciei) const = 0;

        // - Diffusion coefficient of specie ith in the mixture for patch [m^2/ s]
        virtual tmp<scalarField> Dimix( const label speciei, const label patchi) const = 0;

        // - Thermal diffusion coefficient of specie ith in the mixture [m ^2/ s]
        virtual tmp<volScalarField> DimixT(const label speciei) const = 0;
        // - Thermal diffusion coefficient of specie ith in the mixture for patch [m^2/ s]
        virtual tmp<scalarField> DimixT( const label speciei, const label patchi) const = 0;

        // - Binary diffusion coefficient of species ith and species jth in the mixture [m^2/ s] - JH
        virtual tmp<volScalarField> Dij(const label speciei, const label speciej) const = 0;

        //- Enthalpy/Internal energy of specie ith [J/kg]
        virtual tmp<volScalarField> hei(label speciei) const = 0;

        //- Enthalpy/Internal energy of specie ith for patch [J/kg]
        virtual tmp<scalarField> hei(label speciei, const label patchi) const = 0;

        //- Molecular weight of mixture [kg/kmol]
        virtual tmp<volScalarField> Wmix() const = 0;       
        // 

};


/*---------------------------------------------------------------------------*\
                     Class fluidReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidReactionThermo::implementation
:
    virtual public fluidReactionThermo
{

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
