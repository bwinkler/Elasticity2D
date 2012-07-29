example3_setup

dispDisplacement(ds, zExact', 0.1);
eps = 1E-8;
A0 = gf_mesh_fem_get(mfd, 'eval', {1} );
muExact = ds.mVec;

profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLS', 'Adjoint Stiffness');

[muComp, hist, cost, muHist] = is.solve();

profile off;

profsave( profile('info'), mfilename() );

disp( sprintf('\n%d optimization iterations\n ', size(muHist,1)) );

dispMuComparison(ds, muExact, muComp');

printMuComparison(ds, muExact, muComp', mfilename(), [1, 1.4] );

