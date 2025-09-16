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

#include "DTMMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::DTMMultiComponentMixture<ThermoType>::
DTMMultiComponentMixture
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
    linearityMk_(this->Y().size()),
    epsilonOverKbMk_(this->Y().size()),
    sigmaMk_(this->Y().size()),
    miuiMk_(this->Y().size()),
    polarMk_(this->Y().size()),
    ZrotMk_(this->Y().size()),
    CpCoeffTableMk_(this->Y().size()),

    EPSILONijOVERKB_(this->Y().size()),
    DELTAij_(this->Y().size()),
    Mij_(this->Y().size()),
    SIGMAij_(this->Y().size())
    //
{
    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = this->specieThermos()[i].W();
        linearityMk_[i]     = this->specieThermos()[i].linearity();
        epsilonOverKbMk_[i] = this->specieThermos()[i].epsilonOverKb();
        sigmaMk_[i]         = this->specieThermos()[i].sigma();
        miuiMk_[i]          = this->specieThermos()[i].miui();
        polarMk_[i]         = this->specieThermos()[i].polar();
        ZrotMk_[i]          = this->specieThermos()[i].Zrot();        
        CpCoeffTableMk_[i]  = this->specieThermos()[i].CpCoeffTable();
    }

    List<scalar> nEPSILONijOVERKB(this->Y().size());
    List<scalar> nDELTAij(this->Y().size());
    List<scalar> nMij(this->Y().size());
    List<scalar> nSIGMAij(this->Y().size());    

    // Dip = dipole moment 
    const scalar DipMin = 1e-20;

    forAll(EPSILONijOVERKB_, i)
    {
        forAll(nEPSILONijOVERKB, j)
        {
            nMij[j] = 
               1/(1/this->specieThermos()[i].W() + 1/this->specieThermos()[j].W());

            nEPSILONijOVERKB[j] = 
               sqrt(this->specieThermos()[i].epsilonOverKb()*this->specieThermos()[j].epsilonOverKb());

            nSIGMAij[j] = 
               0.5*(this->specieThermos()[i].sigma() + this->specieThermos()[j].sigma());

            // calculate coeficient xi
            scalar xi = 1.0; 
            if ((this->specieThermos()[i].miui() < DipMin) && (this->specieThermos()[j].miui() > DipMin)) 
            {
                // miui_j > DipMin > miui_i --> j is polar, i is nonpolar              
                nDELTAij[j] = 0;
                xi = 1.0 + 
                     this->specieThermos()[i].polar()*pow(this->specieThermos()[j].miui(), 2)* 
                     sqrt(this->specieThermos()[j].epsilonOverKb()/this->specieThermos()[i].epsilonOverKb())/
                     (4*pow(this->specieThermos()[i].sigma(), 3)*
                            (
                             1e+19*this->specieThermos()[j].Kb()*this->specieThermos()[j].epsilonOverKb()*
                             pow(this->specieThermos()[j].sigma(), 3)
                            )
                     );
            }
            else if ((this->specieThermos()[i].miui() > DipMin) && (this->specieThermos()[j].miui() < DipMin))
            {
                // miui_j < DipMin < miui_i --> i is polar, j is nonpolar
                nDELTAij[j] = 0;
                xi = 1.0 +
                     this->specieThermos()[j].polar()*pow(this->specieThermos()[i].miui(), 2)*
                     sqrt(this->specieThermos()[i].epsilonOverKb()/this->specieThermos()[j].epsilonOverKb())/
                     (4*pow(this->specieThermos()[j].sigma(), 3)*   
                            (
                             1e+19*this->specieThermos()[i].Kb()*this->specieThermos()[i].epsilonOverKb()*
                             pow(this->specieThermos()[i].sigma(), 3)
                            )
                     ); 
            }
            else 
            {            
                xi = 1.0;
                nDELTAij[j] =
                   0.5*(this->specieThermos()[i].miui()*this->specieThermos()[j].miui())/
                   (nEPSILONijOVERKB[j]*1e+19*this->specieThermos()[j].Kb()*pow(nSIGMAij[j], 3));
            }
           
            nEPSILONijOVERKB[j] = nEPSILONijOVERKB[j]*pow(xi, 2);
            nSIGMAij[j] = nSIGMAij[j]*pow(xi, -1/6);
        }

        EPSILONijOVERKB_[i] = nEPSILONijOVERKB;
        Mij_[i]             = nMij;
        SIGMAij_[i]         = nSIGMAij;
        DELTAij_[i]         = nDELTAij;
    }
    // End of pre-calculation for kinetic model
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::DTMMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::DTMMultiComponentMixture<ThermoType>::cellThermoMixture
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
    mixture_.updateTRANS
    (
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
    // 

    return mixture_;
}


template<class ThermoType>
const typename
Foam::DTMMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::DTMMultiComponentMixture<ThermoType>::patchFaceThermoMixture
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
    mixture_.updateTRANS
    (
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
    //

    return mixture_;
}






// ************************************************************************* //
