example1_setup
zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

% [zNoise, errNorm] = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );

% zNoise = ds.Q * zNoise + ds.U0';


% Use exact data
pEst = estimateP(ds, zExact );
Z = [zExact; pExact];

A0 = gf_mesh_fem_get(mfd, 'eval', {2.5} );
muExact = ds.mVec;
eps = 1E-9;
%profile on;
is = InverseSolver( ds, A0', Z, eps, 'OLS', 'Hybrid');
[muComp, histOLS, costOLS, muHistOLS] = is.solve();
%profile off;
%profsave( profile('info'), mfilename() );

disp(sprintf('%d optimization iterations\n', size(histOLS,1) ));

dispMuComparison( ds, muExact, muComp');

