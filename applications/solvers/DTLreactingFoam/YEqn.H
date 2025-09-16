tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,Yi_h)")
    )
);

reaction->correct();

forAll(Y, i)
{

    // calculate correction terms
    sumDiffusionCorrections1 *= 0.;
    sumDiffusionCorrections2 *= 0.;
    forAll(Y, k)
    {
        sumDiffusionCorrections1 += fvc::laplacian(Y[i]*rho*Dimix[k], Y[k]);
        sumDiffusionCorrections2 += fvc::laplacian(Y[i]*rho*Dimix[k]*Y[k]/Wmix, Wmix);
    }
    

    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          - fvm::laplacian(rho*Dimix[i], Yi)
         ==
            reaction->R(Yi)
          + fvc::laplacian(rho*Dimix[i]*Yi/Wmix, Wmix) 
          - sumDiffusionCorrections1 
          - sumDiffusionCorrections2     
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}
composition.normalise();

//-Calculate diffusion velocity of inert species
volVectorField YVt(0.0*YVi[0]);

forAll(YVi, i)
{
    //- Update diffusion velocity for all species
    YVi[i] = (-(Dimix[i]/Wmix)*Wmix*fvc::grad(Y[i])-(Dimix[i]/Wmix)*Y[i]*fvc::grad(Wmix));
    if (composition.solve(i))
    {
        YVt += YVi[i];
    }
}
YVi[composition.defaultSpecie()] = -YVt;
