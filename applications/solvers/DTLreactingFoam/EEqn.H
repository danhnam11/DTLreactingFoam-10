{
    volScalarField& he = thermo.he();

    //-Reset and then calculate Heat by diffusion terms-for each time step
    sumHeatDiffusion *= 0.;
    sumHeatDiffusion2 *= 0.;

    forAll(Y, k)
    {
        sumHeatDiffusion += fvc::laplacian((kappa/Cp)*hei[k], Y[k]);
        sumHeatDiffusion2 += fvc::div(hei[k]*rho*YVi[k]);
    }

    heatReactions = reaction->Qdot();

    
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      - fvm::laplacian(kappa/Cp, he)
     ==
        reaction->Qdot()
      - sumHeatDiffusion  
      - sumHeatDiffusion2      
      + fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    #include "calculateAndUpdateCoTHERMFlags.H" 
    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
