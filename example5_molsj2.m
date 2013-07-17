example5_setup

noisePercent = 0.001;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

zNoise = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );
zNoise = ds.Q * zNoise + ds.U0';

pEst = estimateP(ds, zExact );

Z = [zExact; pEst];
%dispDisplacement(ds, zExact', 0.01);

%return 
A0 = gf_mesh_fem_get(mfd, 'eval', {0.2} );

muExact = ds.mVec;
%A0 = muExact;

%smooth
%eps = 1E-5;

% 0.1% noise
%eps = 9E-3;

eps = 1E-9

%BV
%eps = 1E-9;

disp('Starting inverse solver...');
%profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLSJ2', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
%profile off;

profsave( profile('info'), mfilename() );

%disp(sprintf('\n%d optimization iterations\n', size(hist,1) ));

dispMuComparison( ds, muExact, muComp');

%printMuComparison(ds, muExact, muComp', 'example5_ee', [0.5, 1.5] );

% printMuComparison(ds, muExact, muComp', muHist, 'example5_ee', [0.5, 1.5] );

%printMuComparisonStepped(ds, muExact, muComp', muHist, 'example5_molsj2', [0.5 1.5], [0.5 1.5])

%save('example5_molsj1.mat')
