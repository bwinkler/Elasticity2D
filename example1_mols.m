example1_setup

dispDisplacement(ds, zExact', 0.25);

A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );
muExact = ds.mVec;

%smooth
%eps = 5E-8;

%BV
eps = 1E-9;

disp('Starting inverse solver...');
tic
is = InverseSolver( ds, A0', Z, eps, 'MOLS', 'Adjoint Stiffness');
[muComp, hist, cost, muHist] = is.solve();
toc

disp(sprintf('%d optimization iterations\n', size(hist,1) ));

dispMuComparison( ds, muExact, muComp');

%print -zbuffer -r600 -dtiff example2.tiff

