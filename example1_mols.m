example1_setup

noisePercent = 0.01;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

[zNoise, errNorm] = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );

zNoise = ds.Q * zNoise + ds.U0';


% Use exact data
pEst = estimateP(ds, zExact );
Z = [zExact; pEst];
eps = 1E-8;

% Use noise data directly
% pEst = estimateP(ds, zNoise );
% Z = [zNoise; pEst];
%eps = 2E-4;

% Use smoothed data
% res = @(a)  errNorm^2 - norm( zNoise - smoothData1(zNoise, a, ds))^2;
% alpha = secant(res, 2.6E-2, 3.7E-3, 1E-15, 1E-15, 100000000)
% zSmooth = smoothData1( zNoise, alpha, ds );
% pEst = estimateP(ds, zSmooth);
% Z = [zSmooth; pEst];
%norm( zExact - zSmooth)
%dispDisplacementComparison( ds, zNoise', zSmooth', 1);
%eps = 2E-3;

%return 
A0 = gf_mesh_fem_get(mfd, 'eval', {1} );
muExact = ds.mVec;

disp('Starting inverse solver...');
%profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
%profile off;

%profsave( profile('info'), mfilename() );

disp(sprintf('\n%d optimization iterations\n', size(hist,1) ));

dispMuComparison( ds, muExact, muComp');

%printMuComparison(ds, muExact, muComp', mfilename(), [2.2, 2.8] );

