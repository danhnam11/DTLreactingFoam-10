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

\*---------------------------------------------------------------------------*/

#include "FTMMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FTMMultiComponentMixture<ThermoType>::
FTMMultiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    mixture_("mixture", this->specieThermos()[0]),

    // for mixture DTM - Nam
    ListW_(this->Y().size()),
    muCoeffsMk_(this->Y().size()),
    kappaCoeffsMk_(this->Y().size()),
    DijCoeffsMk_(this->Y().size())
    //
{
    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]         = this->specieThermos()[i].W();
        muCoeffsMk_[i]    = this->specieThermos()[i].coeffs().muCoeffs();
        kappaCoeffsMk_[i] = this->specieThermos()[i].coeffs().kappaCoeffs();
        DijCoeffsMk_[i]   = this->specieThermos()[i].coeffs().DijCoeffs();
    }
    // End of pre-calculation for kinetic model

    // update coeffs for mixture properties here to save CPU time
    // since these coffs are constant
    mixture_.updateTRANSFitCoeff
    (
        muCoeffsMk_, kappaCoeffsMk_, DijCoeffsMk_
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::FTMMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::FTMMultiComponentMixture<ThermoType>::cellThermoMixture
(
    const label celli
) const
{
    mixture_ = this->Y()[0][celli]*this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ += this->Y()[i][celli]*this->specieThermos()[i];
    }

    //- Update coefficients for mixture DTM - Nam
    //- List of secies mole and mass fraction 
    List<scalar> X(this->Y().size()); 
    List<scalar> Y(this->Y().size()); 
    scalar sumXb = 0.0;  
    forAll(X, i)
    {
        sumXb = sumXb + this->Y()[i][celli]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(X, i)
    {
        X[i] = (this->Y()[i][celli]/ListW_[i])/sumXb;
        Y[i] = this->Y()[i][celli];
        if(X[i] <= 0) { X[i] = 0; }
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }
    
    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    // Update coefficients for mixture DTM - Nam
    mixture_.updateTRANS(Y, X, ListW_);
    // 

    return mixture_;
}


template<class ThermoType>
const typename
Foam::FTMMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::FTMMultiComponentMixture<ThermoType>::patchFaceThermoMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ =
        this->Y()[0].boundaryField()[patchi][facei]
       *this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ +=
            this->Y()[i].boundaryField()[patchi][facei]
           *this->specieThermos()[i];
    }

    //- Update coefficients for mixture - DTM
    //- List of secies mole and mass fraction 
    List<scalar> X(this->Y().size());
    List<scalar> Y(this->Y().size());
    scalar sumXb = 0.0;
    forAll(X, i)
    {
        sumXb = sumXb + this->Y()[i].boundaryField()[patchi][facei]/ListW_[i];
    }
    if (sumXb == 0){ sumXb = 1e-30;}

    forAll(X, i)
    {
        X[i] = (this->Y()[i].boundaryField()[patchi][facei]/ListW_[i])/sumXb;
        Y[i] = this->Y()[i].boundaryField()[patchi][facei];
        if(X[i] <= 0) { X[i] = 0; } 
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }

    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    // Update coefficients for mixture DTM - Nam
    mixture_.updateTRANS(Y, X, ListW_);
    //

    return mixture_;
}



// ************************************************************************* //
