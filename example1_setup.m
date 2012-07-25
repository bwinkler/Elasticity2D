clear all;
gf_workspace('clear all');
warning off;

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

% For totally incompressible case.
%ld = 'Incompressible';

f1 = '1+0.1*x';
f2 = '1+0.1*y';


dirBound1 = '0.1*x';
dirBound2 = '0.1*y';

neuBound1 = '1';
neuBound2 = '1';


disp('Starting direct solver...');
tic
ds = DirectSolver( m, mu, ld, f1, f2, mfu, mfd, mfp, mim, ...
                   TOP, dirBound1, dirBound2, ...
                   REST, neuBound1, neuBound2,...
                   'H1Semi');
Z = ds.solve();
toc;

zExact = Z(1:ds.udof);
pExact = Z((ds.udof+1):length(Z));

