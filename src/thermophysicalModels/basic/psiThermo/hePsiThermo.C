/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "hePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// calculate thermo variables using original transport models in OF
template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] =
        transportMixture.kappa(pCells[celli], TCells[celli]);

    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
    }
}


// calculate thermo variables using DTM
// only use for constructor since it is relatively expensive
template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::initialize()
{
    Info << "[Expensive!] using Detailed Transport Model WITHOUT coTHERM " << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = Wmix_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();

        forAll(Dimix_, i)
        {
            Dimix_[i].primitiveFieldRef()[celli] 
          = transportMixture.Dimix(i, pCells[celli], TCells[celli]);

            DimixT_[i].primitiveFieldRef()[celli]
          = transportMixture.DimixT(i, pCells[celli], TCells[celli]);

            heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].HE(pCells[celli], TCells[celli]);
        }
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = Wmix_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();

                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
                pWmix[facei]  = thermoMixture.W();

                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);                  
                }
            }
        }
    }
}


// calculate thermo variables using DTM with coTHERM algorithm 
// This can significantly reduce computational time
template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateUsingCoTHERM()
{
    Info << "[Note!] using Detailed Transport Models + CoTHERM" << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = Wmix_.primitiveFieldRef(); 

    // For CoTHERMCount
    scalarField& CountCells = this->coTHERMStepCount_.primitiveFieldRef();

    scalarField& TCellsOld = this->T_.oldTime().primitiveFieldRef();
    scalarField& PCellsOld = this->p_.oldTime().primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        this->residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli];
        scalar coDeltaT = mag(TCells[celli] - TCellsOld[celli]);
        scalar coDeltaP = mag(pCells[celli] - PCellsOld[celli]);

        const bool exceedCount = (CountCells[celli] >= this->maxCoTHERMStepCount_);

        if 
        ( 
            (coDeltaT <= this->epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            (coDeltaP <= this->epsilonP_) &&
            (!exceedCount)              
        )
        {
            // If the temperature is unchanged, don't need to calculate 
            // thermophysical properties, just copy from old time fields 
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 1.0;
            // start counter     
            CountCells[celli] += 1.0;

            psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
            WmixCells[celli]  = thermoMixture.W();

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = Dimix_[i].oldTime()[celli];
    
                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];
    
                heList_[i].primitiveFieldRef()[celli]
              = heList_[i].oldTime()[celli];
            }
        }
        else if 
        (
            (coDeltaT <= this->epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            !(coDeltaP <= this->epsilonP_) &&
            (!exceedCount)               
        )
        {
            // only re-calculate Dimix
            // for other thermophysical properties, copy from old time
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 0.5;
            // reset counter
            CountCells[celli] = 0.0; 

            psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
            WmixCells[celli]  = thermoMixture.W();

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
    
                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];
    
                heList_[i].primitiveFieldRef()[celli]
              = heList_[i].oldTime()[celli];
            }
        }
        else
        {
            // otherwise, calculate all thermophysical properties again 
            this->coTHERMstatus_.primitiveFieldRef()[celli] = 0.0;
            // reset counter
            CountCells[celli] = 0.0; 

            CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
            CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
            psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
            WmixCells[celli]  = thermoMixture.W();  
            muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
            kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
    
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
    
                DimixT_[i].primitiveFieldRef()[celli]
              = transportMixture.DimixT(i, pCells[celli], TCells[celli]);
    
                heList_[i].primitiveFieldRef()[celli]
              = specieThermos_[i].HE(pCells[celli], TCells[celli]);
            }
        }
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = Wmix_.boundaryFieldRef();
    volScalarField::Boundary& CountBf  = this->coTHERMStepCount_.boundaryFieldRef();

    // old fields
    volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& pBfOld = this->p_.oldTime().boundaryFieldRef();   
    volScalarField::Boundary& CpBfOld = this->Cp_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& CvBfOld = this->Cv_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& muBfOld = this->mu_.oldTime().boundaryFieldRef();    
    volScalarField::Boundary& kappaBfOld = this->kappa_.oldTime().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];
        fvPatchScalarField& pCount  = CountBf[patchi];

        // old fields
        fvPatchScalarField& pTOld = TBfOld[patchi]; 
        fvPatchScalarField& ppOld = pBfOld[patchi];         
        fvPatchScalarField& pCpOld = CpBfOld[patchi];
        fvPatchScalarField& pCvOld = CvBfOld[patchi];        
        fvPatchScalarField& pmuOld = muBfOld[patchi];        
        fvPatchScalarField& pkappaOld = kappaBfOld[patchi];
        //

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);

                const bool exceedCountP = (pCount[facei] >= this->maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    (coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                                           
                )
                {
                    // don't need to re-calculate thermophysical properties
                    // just copy from old time fields
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0; 
                    // start counter
                    pCount[facei] += 1.0;

                    phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                    pWmix[facei]  = thermoMixture.W();                    
    
                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else if 
                (
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    !(coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                         
                )
                {
                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.5;
                    // reset counter
                    pCount[facei] = 0.0;

                    phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                    pWmix[facei]  = thermoMixture.W();                    
    
                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate all thermophysical properties again 
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;
                    // reset counter
                    pCount[facei] = 0.0;

                    phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                    pWmix[facei]  = thermoMixture.W();                    
    
                    pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                    pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                    pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                    pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.DimixT(i, pp[facei], pT[facei]);
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);
                    }
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);

                const bool exceedCountP = (pCount[facei] >= this->maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    (coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                                          
                )
                {
                    // don't need to re-calculate thermophysical properties
                    // just copy from old time fields
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 1.0;
                    // start counter
                    pCount[facei] += 1.0;

                    pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]); 
                    pWmix[facei]  = thermoMixture.W();                                                    

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else if 
                (
                    (coDeltaTp <= this->epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    !(coDeltaPp <= this->epsilonP_) &&
                    (!exceedCountP)                         
                )
                {
                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.5;
                    // reset counter
                    pCount[facei] = 0.0;

                    pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]); 
                    pWmix[facei]  = thermoMixture.W();                                                    

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate all thermophysical properties again     
                    this->coTHERMstatus_.boundaryFieldRef()[patchi][facei] = 0.0;
                    // reset counter
                    pCount[facei] = 0.0;
                            
                    pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]); 
                    pWmix[facei]  = thermoMixture.W();                                                    

                    pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                    pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                    pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                    pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dimix(i, pp[facei], pT[facei]);
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = transportMixture.DimixT(i, pp[facei], pT[facei]);
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);                  
                    }
                }
            }
        }
    }
}


// calculate thermo variables using DTM with coTHERM algorithm, 
// but only using T criterion, only internally in our Lab
// This also can significantly reduce computational time
template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateUsingCoTHERMOnlyT()
{
    Info << "[Note!] using Detailed Transport Model + CoTHERM only T " << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    scalarField& WmixCells = Wmix_.primitiveFieldRef(); //

    scalarField& TCellsOld = this->T_.oldTime().primitiveFieldRef(); //

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        this->residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli]; //

        if ( mag(TCells[celli] - TCellsOld[celli] ) < this->epsilonOnlyT_ )   
        {
            // If the temperature is unchanged, don't need to calculate 
            // thermophysical properties, just copy from old time fields 
            psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
            WmixCells[celli]  = thermoMixture.W();

            CpCells[celli] = this->Cp_.oldTime()[celli];
            CvCells[celli] = this->Cv_.oldTime()[celli];
            muCells[celli] = this->mu_.oldTime()[celli];
            kappaCells[celli] = this->kappa_.oldTime()[celli];
    
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = Dimix_[i].oldTime()[celli];
    
                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];
    
                heList_[i].primitiveFieldRef()[celli]
              = heList_[i].oldTime()[celli];
            }
        }
        else
        {
            // otherwise, calculate thermophysical properties again 
			CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
			CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
			psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
            WmixCells[celli]  = thermoMixture.W(); 	
			muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
			kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
	
			forAll(Dimix_, i)
			{
				Dimix_[i].primitiveFieldRef()[celli] 
			  = transportMixture.Dimix(i, pCells[celli], TCells[celli]);
	
				DimixT_[i].primitiveFieldRef()[celli]
			  = transportMixture.DimixT(i, pCells[celli], TCells[celli]);
	
				heList_[i].primitiveFieldRef()[celli]
			  = specieThermos_[i].HE(pCells[celli], TCells[celli]);
			}
        }
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = Wmix_.boundaryFieldRef(); // 

	// old fields
    volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& CpBfOld = this->Cp_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& CvBfOld = this->Cv_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& muBfOld = this->mu_.oldTime().boundaryFieldRef();    
    volScalarField::Boundary& kappaBfOld = this->kappa_.oldTime().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        fvPatchScalarField& pWmix = WmixBf[patchi]; // 

        // old fields
        fvPatchScalarField& pTOld = TBfOld[patchi]; 
        fvPatchScalarField& pCpOld = CpBfOld[patchi];
        fvPatchScalarField& pCvOld = CvBfOld[patchi];        
        fvPatchScalarField& pmuOld = muBfOld[patchi];        
        fvPatchScalarField& pkappaOld = kappaBfOld[patchi];
		//

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei]; //

                if ( mag(pT[facei] - pTOld[facei] ) < this->epsilonOnlyT_ )
                {
                    // If the temperature is unchanged, don't need to calculate 
                    // thermophysical properties, just copy from old time fields 
                    phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                    pWmix[facei]  = thermoMixture.W();                    
    
                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate thermophysical properties again                     
					phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                    pWmix[facei]  = thermoMixture.W();                    
	
					pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
					pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
					pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
					pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
	
					forAll(Dimix_, i)
					{
						Dimix_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.Dimix(i, pp[facei], pT[facei]);
	
						DimixT_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.DimixT(i, pp[facei], pT[facei]);
	
						heList_[i].boundaryFieldRef()[patchi][facei]
					  = specieThermos_[i].HE(pp[facei], pT[facei]);
					}
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                this->residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei]; //
				
                if ( mag(pT[facei] - pTOld[facei] ) < this->epsilonOnlyT_ )
                {
                    // If the temperature is unchanged, don't need to calculate 
                    // thermophysical properties, just copy from old time fields 
                    pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]); 
                    pWmix[facei]  = thermoMixture.W();                                                    

                    pCp[facei] = pCpOld[facei];
                    pCv[facei] = pCvOld[facei];
                    pmu[facei] = pmuOld[facei];
                    pkappa[facei] = pkappaOld[facei];
    
                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];
    
                        heList_[i].boundaryFieldRef()[patchi][facei]
                      = heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate thermophysical properties again     
					pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                    ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]); 
                    pWmix[facei]  = thermoMixture.W();                                                    

					pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
					pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
					pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
					pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
	
					forAll(Dimix_, i)
					{
						Dimix_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.Dimix(i, pp[facei], pT[facei]);
	
						DimixT_[i].boundaryFieldRef()[patchi][facei]
					  = transportMixture.DimixT(i, pp[facei], pT[facei]);
	
						heList_[i].boundaryFieldRef()[patchi][facei]
					  = specieThermos_[i].HE(pp[facei], pT[facei]);                  
					}
                }
            }
        }
    }
}



// calculate thermo variables using DTM + generate binary diffusivity coeffs
// only use this function in preprocessing steps for FTM model
template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateTransportPreProcessing()
{
    Info << "[Expensive!] using Detailed Transport Model for preprocessing FTM " << endl;    
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    scalarField& WmixCells = Wmix_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] = transportMixture.kappa(pCells[celli], TCells[celli]);
        WmixCells[celli]  = thermoMixture.W();

        forAll(Dimix_, i)
        {
            Dimix_[i].primitiveFieldRef()[celli] 
          = transportMixture.Dimix(i, pCells[celli], TCells[celli]);

            DimixT_[i].primitiveFieldRef()[celli]
          = transportMixture.DimixT(i, pCells[celli], TCells[celli]);

            heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].HE(pCells[celli], TCells[celli]);
        }

        // for binary diffusion coefficients
        forAll(Dij_, i)
        {
            forAll(Dij_[i],j)
            {
                Dij_[i][j].primitiveFieldRef()[celli]
              = transportMixture.Dij(i, j, pCells[celli], TCells[celli]);
            }
        }
     
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = Wmix_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);

                pWmix[facei]  = thermoMixture.W();
                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }

                // for binary diffusion coefficients
                forAll(Dij_, i)
                {
                    forAll(Dij_[i], j)
                    {
                        Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dij(i, j, pp[facei], pT[facei]);
                    }
                }

            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);

                pWmix[facei]  = thermoMixture.W();                
                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = transportMixture.DimixT(i, pp[facei], pT[facei]);

                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);                  
                }

                // for binary diffusion coefficients
                forAll(Dij_, i)
                {
                    forAll(Dij_[i], j)
                    {
                        Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = transportMixture.Dij(i, j, pp[facei], pT[facei]);
                    }
                }
                
            }
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::hePsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    // Nam 
    Dimix_(MixtureType::specieThermos().size()),  
    DimixT_(MixtureType::specieThermos().size()), 
    Dij_(MixtureType::specieThermos().size()),    
    heList_(MixtureType::specieThermos().size()), 
    Wmix_
    (
        IOobject
        (
            "Wmix",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Wmix", dimensionSet(1, 0, 0, 0, -1), 0.)
    ),
    specieThermos_(this->specieThermos()),
    mode_(0) //    
{

    // for diffusivity models in DTM - Nam
    forAll(Dimix_, i)
    {
        Dimix_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    this->phasePropertyName("thermo:Dimix"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionSet(0, 2, -1, 0, 0)
            )
        );
    }

    forAll(DimixT_, i)
    {
        DimixT_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    this->phasePropertyName("thermo:DimixT"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionSet(1, -1, -1, 0, 0)
            )
        );
    }
    forAll(Dij_, i)
    {
        Dij_[i] = PtrList<volScalarField>(MixtureType::specieThermos().size());

        forAll(Dij_[i], j)
        {
            Dij_[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        this->phasePropertyName("thermo:Dij"),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionSet(0, 2, -1, 0, 0)
                )
            );
        }
    }

    forAll(heList_, i)
    {
        heList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hei",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                this->he_.dimensions()
            ) 
        );
    }

    // check and setup coTHERM mode - Nam
    if 
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "DetailedModels";        
        mode_ = 1;
    }
    else if 
    (
        this->usingDetailedTransportModel_ && 
        this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "coTHERM";        
        mode_ = 2;
    }
    else if
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        !this->usingCoTHERMOnlyT_ &&
        this->usingPreProcessingFTM_        
    )
    {
        //modeName_ = "preprocessingCoTHERM";
        mode_ = 3;
    }
    else if 
    (
        this->usingDetailedTransportModel_ && 
        !this->usingCoTHERM_ && 
        this->usingCoTHERMOnlyT_ && 
        !this->usingPreProcessingFTM_       
    )
    {
        //modeName_ = "CoTHERMonlyT";
        mode_ = 4;
    }
    else
    {
        //modeName_ = "originalOpenFOAM";
        mode_ = 0;
    }
    // Nam 

    // calculate();  // original 
    Info << "[This is in constructor] thermophysical properties initialization! "  << endl;
    initialize(); // Nam - for DTM
    Info << "[This is in constructor] end of thermophysical properties initialization! "  << endl;

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::~hePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    //calculate(); // original

    // Nam
    switch(mode_)
    {
        case 1 : 
            initialize();            
            break; 
            
        case 2 : 
            calculateUsingCoTHERM();
            break;      
            
        case 3 : 
            calculateTransportPreProcessing();
            break; 

        case 4 : 
            calculateUsingCoTHERMOnlyT();
            break; 
            
        default:
            calculate();
            break; 
    }
    // Nam 

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
