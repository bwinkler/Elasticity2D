example2_setup

eps = 1E-6;
A0 = gf_mesh_fem_get(mfd, 'eval', {0.15} );
muExact = ds.mVec;

%profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLSJ2', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
%profile off;


%profsave( profile('info'), mfilename() );

disp( sprintf('\n%d optimization iterations\n ', size(muHist,1)) );

%printMuComparison(ds, muExact, muComp', mfilename(), [0.2, 0.6] );

dispMuComparison(ds, muExact, muComp');

%save('example2_molsj1.mat')

%printMuComparisonStepped(ds, muExact, muComp', muHist, 'example2_molsj1', [0.2 0.6], [0.1.5 1])

