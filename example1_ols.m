example1_setup

A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );
muExact = ds.mVec;
eps = 1E-6;
tic
is = InverseSolver( ds, A0', Z, eps, 'OLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
toc

dispMuComparison( ds, muExact, muComp');

