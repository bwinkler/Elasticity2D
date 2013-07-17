example3_setup

noisePercent = 0.01;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

[zNoise, noiseError] = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );

zNoise = ds.Q * zNoise + ds.U0';

% use exact solution for Z
pEst = estimateP( ds, zExact);
Z = [ zExact; pEst ];
eps = 1E-6;

% Use noise data directly
% pEst = estimateP(ds, zNoise );
% Z = [zNoise; pEst];
% eps = 1E-4;

% use smoothed data
% res = @(a)  noiseError^2 - norm( zNoise - smoothData2(zNoise, a, ds))^2;
% alpha = secant(res, 2.6E-2, 3.7E-3, 1E-15, 1E-15, 100000)
% zSmooth = smoothData2( zNoise, alpha, ds );
% pEst = estimateP(ds, zSmooth);
% Z = [zSmooth; pEst];
% eps = 1E-4;
% smoothError = norm( zExact - zSmooth);
% smoothEff = noiseError - smoothError
% dispDisplacementComparison( ds, zNoise', zSmooth', 0.1);

%dispDisplacement(ds, zExact', 0.1);

A0 = gf_mesh_fem_get(mfd, 'eval', {0.5} );
%A0 = A0 + rand(size(A0));
muExact = ds.mVec;

%profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLSJ2', 'Adjoint Stiffness');

[muComp, hist, cost, muHist] = is.solve();

%profile off;

%profsave( profile('info'), mfilename() );

disp( sprintf('\n%d optimization iterations\n ', size(muHist,1)) );

dispMuComparison(ds, muExact, muComp', [15 30]);

%printMuComparisonStepped(ds, muExact, muComp', muHist, 'example3_molsj2', [0.5 1.3], [0.5 1.3], [1,4,8],[15 30])
%save('example3_molsj2.mat')
%printMuComparison(ds, muExact, muComp', mfilename(), [1, 1.4] );

