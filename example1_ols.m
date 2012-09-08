example1_setup

A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );
muExact = ds.mVec;
eps = 1E-8;
tic
is = InverseSolver( ds, A0', Zt, eps, 'OLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
toc

disp(sprintf('%d optimization iterations\n', size(hist,1) ));

dispMuComparison( ds, muExact, muComp');

