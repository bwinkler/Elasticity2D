example4_setup

noisePercent = 0.001;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

zNoise = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );
zNoise = ds.Q * zNoise + ds.U0';

pEst = estimateP(ds, zNoise );

Z = [zNoise; pEst];
%dispDisplacement(ds, zExact', 0.01);

%return 
A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );

muExact = ds.mVec;
%A0 = muExact;

%smooth
%eps = 1E-5;

% 0.1% noise
eps = 9E-3;

%BV
%eps = 1E-9;

disp('Starting inverse solver...');
%profile on;
is = InverseSolver( ds, A0', Z, eps, 'EE', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
%profile off;

%profsave( profile('info'), mfilename() );

%disp(sprintf('\n%d optimization iterations\n', size(hist,1) ));

dispMuComparison( ds, muExact, muComp');

%printMuComparison(ds, muExact, muComp', mfilename(), [2.2, 2.8] );

