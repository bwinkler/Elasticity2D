
N = 1;
trialsMOLS = zeros(1,N);
for i = [1:N]
example1_setup

%noisePercent = 0.01;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

%[zNoise, errNorm] = addNoise( ds.Q' * zExact, noisePercent, 'uniform' );

%zNoise = ds.Q * zNoise + ds.U0';


% Use exact data
% Z = [zExact; pExact];
% eps = 1E-8;

% Use estimate
pEst = estimateP(ds, zExact );
Z = [zExact; pEst];
eps = 1E-6;

% Use interpolated
% pEst = estimateP(ds, zInterp );
% Z = [zInterp; pEst];
% eps = 1E-6;

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
%dispDisplacementComparison( ds, zNoise', zSmooth', 1);
%eps = 2E-3;

% norm( zExact - zInterp)
% dispDisplacementComparison( ds, zExact', zInterp', 1);
% return 

A0 = gf_mesh_fem_get(mfd, 'eval', {1} );
%A0 = A0 + rand(size(A0));
muExact = ds.mVec;

disp('Starting inverse solver...');
%profile on;
is = InverseSolver( ds, A0', Z, eps, 'MOLSJ1', 'Adjoint');
[muComp, histMOLS, costMOLS, muHistMOLS] = is.solve();
%profile off;

%profsave( profile('info'), mfilename() );

disp(sprintf('\n%d optimization iterations\n', size(histMOLS,1) ));

%trialsMOLS(i) = norm(muExact - muComp');
end

%mean(trialsMOLS)
%save('trialsMOLS.mat', 'trialsMOLS')
dispMuComparison( ds, muExact, muComp');

%printMuComparison(ds, muExact, muComp', mfilename(), [2.2, 2.8] );

%printMuComparison(ds, muExact, muComp', muHistMOLS, [2,5,8], mfilename(), [2.2,2.8]);

%printMuComparisonCombined(ds, muExact, muComp', mfilename(), [2.2, 2.8] );

%printMuComparisonStepped(ds, muExact, muComp', muHist, 'example1_molsj1', [2.2 2.8], [1 2.8])
%save('example1_molsj1.mat')

