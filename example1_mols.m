example1_setup

%dispDisplacement(ds, zExact', 0.01);

%return 
A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );
muExact = ds.mVec;

%smooth
eps = 1E-8;

%BV
%eps = 1E-9;

disp('Starting inverse solver...');
profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
profile off;

profsave( profile('info'), mfilename() );

disp(sprintf('\n%d optimization iterations\n', size(hist,1) ));

%dispMuComparison( ds, muExact, muComp');

printMuComparison(ds, muExact, muComp', mfilename(), [2.2, 2.8] );

