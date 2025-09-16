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

#include "heheuPsiThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& heuCells = this->heu_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& TuCells = this->Tu_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    scalarField& WmixCells = Wmix_.primitiveFieldRef(); //  

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

        TuCells[celli] = this->cellReactants(celli).THE
        (
            heuCells[celli],
            pCells[celli],
            TuCells[celli]
        );

        // Nam
        WmixCells[celli]  = thermoMixture.W(); 
        forAll(Dimix_, i)
        {
            Dimix_[i].primitiveFieldRef()[celli] 
          //= transportMixture.Dimix(i, pCells[celli], TCells[celli]);
            = 1; 

            DimixT_[i].primitiveFieldRef()[celli]
          //= transportMixture.DimixT(i, pCells[celli], TCells[celli]);
            = 1;

            heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].HE(pCells[celli], TCells[celli]);
        }
        // Nam 

    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& TuBf = this->Tu_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& heuBf = this->heu().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = Wmix_.boundaryFieldRef(); //  

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pTu = TuBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pheu = heuBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi]; //      

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

                // Nam 
                pWmix[facei]  = thermoMixture.W();                      
                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  //= transportMixture.Dimix(i, pp[facei], pT[facei]);
                    = 1;

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  //= transportMixture.DimixT(i, pp[facei], pT[facei]);
                    = 1; 
 
                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }
                // 

            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

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

                pTu[facei] =
                    this->patchFaceReactants(patchi, facei)
                   .THE(pheu[facei], pp[facei], pTu[facei]);


                // Nam 
                pWmix[facei]  = thermoMixture.W();                      
                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  //= transportMixture.Dimix(i, pp[facei], pT[facei]);
                    = 1;

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  //= transportMixture.DimixT(i, pp[facei], pT[facei]);
                    = 1; 

                    heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }
                // 
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heheuPsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    Tu_
    (
        IOobject
        (
            "Tu",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    heu_
    (
        IOobject
        (
            MixtureType::thermoType::heName() + 'u',
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->volScalarFieldProperty
        (
            MixtureType::thermoType::heName() + 'u',
            dimEnergy/dimMass,
            &MixtureType::cellReactants,
            &MixtureType::patchFaceReactants,
            &MixtureType::thermoMixtureType::HE,
            this->p_,
            this->Tu_
        ),
        this->heuBoundaryTypes()
    ),
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
    specieThermos_(this->specieThermos())   
    // Nam 
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
    // - Nam 

    this->heuBoundaryCorrection(this->heu_);

    calculate();

    this->psi_.oldTime(); // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::~heheuPsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& Tu,
    const labelList& cells
) const
{
    return this->cellSetProperty
    (
        &MixtureType::cellReactants,
        &MixtureType::thermoMixtureType::HE,
        cells,
        UIndirectList<scalar>(this->p_, cells),
        Tu
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& Tu,
    const label patchi
) const
{
    return this->patchFieldProperty
    (
        &MixtureType::patchFaceReactants,
        &MixtureType::thermoMixtureType::HE,
        patchi,
        this->p_.boundaryField()[patchi],
        Tu
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::Tb() const
{
    return this->volScalarFieldProperty
    (
        "Tb",
        dimTemperature,
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::thermoMixtureType::THE,
        this->he_,
        this->p_,
        this->T_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psiu() const
{
    return this->volScalarFieldProperty
    (
        "psiu",
        this->psi_.dimensions(),
        &MixtureType::cellReactants,
        &MixtureType::patchFaceReactants,
        &MixtureType::thermoMixtureType::psi,
        this->p_,
        this->Tu_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psib() const
{
    const volScalarField Tb(this->Tb());

    return this->volScalarFieldProperty
    (
        "psib",
        this->psi_.dimensions(),
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::thermoMixtureType::psi,
        this->p_,
        Tb
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::muu() const
{
    return this->volScalarFieldProperty
    (
        "muu",
        dimDynamicViscosity,
        &MixtureType::cellReactants,
        &MixtureType::patchFaceReactants,
        &MixtureType::transportMixtureType::mu,
        this->p_,
        this->Tu_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::mub() const
{
    const volScalarField Tb(this->Tb());

    return this->volScalarFieldProperty
    (
        "mub",
        dimDynamicViscosity,
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::transportMixtureType::mu,
        this->p_,
        Tb
    );
}


// ************************************************************************* //
