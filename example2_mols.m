example2_setup

eps = 1E-6;
A0 = gf_mesh_fem_get(mfd, 'eval', {0.15} );
muExact = ds.mVec;

tic
is = InverseSolver( ds, A0', Z, eps, 'MOLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
toc

disp( sprintf('\n%d optimization iterations\n ', size(muHist,1)) );

dispMuComparison(ds, muExact, muComp');

%print -zbuffer -r600 -dtiff example2.tiff

