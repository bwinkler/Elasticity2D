clear all;
gf_workspace('clear all');

% Create the mesh
n = 10;
h = 1/(n+1);
m = gf_mesh('cartesian', [0:h:1], [0:h:1] );
%gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% Assign boundary regions
TOP = 1;
REST = 2;
BORDER = 3;
border = gf_mesh_get(m, 'outer faces');

pts = gf_mesh_get(m,'pts');
pidtop = find( abs(pts(2,:) - 1 ) < 1/(2*(n+1)));
top = gf_mesh_get(m, 'faces_from_pid', pidtop);

rest = setdiff(border', top', 'rows')';

gf_mesh_set(m, 'boundary', TOP, top);
gf_mesh_set(m, 'boundary', REST, rest);
gf_mesh_set(m, 'boundary', BORDER, border);
%gf_plot_mesh(m, 'regions', [TOP]);

mfu = gf_mesh_fem(m,2);
gf_mesh_fem_set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));

mfd = gf_mesh_fem(m);
gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_QK(2,1)'));

mfp = gf_mesh_fem(m);
gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_QK(2,0)'));

%mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,2)'));
mim = gf_mesh_im(m, gf_integ('IM_QUAD(2)'));

mu = '2.5+0.25*sin(2*pi*x)';
ld = 1E6;
%ld = 'Incompressible';
%mu = '1.5 + x.*y';
%mu = '1.5+0.5*sin(pi*y)+0.5*sin(pi*x)';
%mu = '2+x.^2+y.^2';
fx = '1+0.1*x';
fy = '1+0.1*y';


dirBoundx = '0';
dirBoundy = '0';

neuBoundx = '1';
neuBoundy = '1';

eps = 5E-8;

tic
ds = DirectSolver( m, mu, ld, fx, fy, mfu, mfd, mfp, mim, ...
                   TOP, dirBoundx, dirBoundy, ...
                   REST, neuBoundx, neuBoundy,...
                   'H1Semi');
Ut = ds.solve();
toc;

u1 = Ut(1:ds.udof);
p = Ut((ds.udof+1):length(Ut));

Z = Ut;
A0 = gf_mesh_fem_get(mfd, 'eval', {2.2} );
mVec = ds.mVec;
is = InverseSolver( ds, A0', Z, eps);
[MUc, hist,cost, MUHist] = is.solve();
toc

Ut = ds.solve();
u = Ut(1:ds.udof);
p = Ut((ds.udof+1):length(Ut));

subplot(1,2,1);
gf_plot( mfd, mVec, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
%colorbar;
title('Exact \mu');
subplot(1,2,2);
gf_plot( mfd, MUc',  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
%colorbar;
title('Computed \mu');


