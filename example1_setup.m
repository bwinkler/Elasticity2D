%clear all;
gf_workspace('clear all');
%warning off;

% Create the mesh
n = 100; 
h = 1/(n+1);


% n0 = 500;
% h0 = 1/(n0+1);
% m0 = gf_mesh('cartesian', [0:h0:1], [0:h0:1]);
% mfu0 =gf_mesh_fem(m0,2);
% gf_mesh_fem_set(mfu0, 'fem', gf_fem('FEM_QK(2,1)'));


m = gf_mesh('cartesian', [0:h:1], [0:h:1] );
%m = gf_mesh('triangles grid', [0:h:1], [0:h:1] );

pts = gf_mesh_get(m,'pts');
pidleft = find( abs(pts(1,:) ) < 1/(2*(n+1)) );
pidright = find( abs(pts(1,:) - 1 ) < 1/(2*(n+1)) );
pidbottom = find( abs(pts(2,:) ) < 1/(2*(n+1)) );
pidtop = find( abs(pts(2,:) - 1 ) < 1/(2*(n+1)) );

% dirBoundFaces = gf_mesh_get(m, 'faces_from_pid', union(pidleft, pidright) );
% neuBoundFaces = gf_mesh_get(m, 'faces_from_pid',...
%                   union( pidtop, pidbottom) );

dirBoundFaces = gf_mesh_get(m, 'faces_from_pid', pidtop );
neuBoundFaces = gf_mesh_get(m, 'faces_from_pid',...
                  union( pidleft, union(pidbottom, pidright)) );
dirBoundID = 1;
neuBoundID = 2;

gf_mesh_set(m, 'boundary', dirBoundID, dirBoundFaces);
gf_mesh_set(m, 'boundary', neuBoundID, neuBoundFaces);

mfu = gf_mesh_fem(m,2);
gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_QK(2,1)'));

mfd = gf_mesh_fem(m);
gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_QK(2,1)'));

mfp = gf_mesh_fem(m);
gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_QK(2,0)'));

% mfu = gf_mesh_fem(m,2);
% gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_PK(2,2)'));

% mfd = gf_mesh_fem(m);
% gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_PK(2,1)'));

% mfp = gf_mesh_fem(m);
% gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_PK(2,0)'));
% %mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,1)'));
% mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(2)'));
mim = gf_mesh_im(m, gf_integ('IM_QUAD(2)'));

mu = '2.5+0.25*sin(2*pi*x)';
ld = 1E6

% For totally incompressible case.
%ld = 'Incompressible';

f1 = '2.3+0.1*x';
f2 = '2.3+0.1*y';

dirBound1 = '0.01*x';
dirBound2 = '0.01*y.^2';

neuBound1 = '0.5+x.^2';
neuBound2 = '0.5+y.^2';

% f1 = '0.001';
% f2 = '0.001';

% dirBound1 = '0';
% dirBound2 = '0';

% neuBound1 = '0.1';
% neuBound2 = '0';


disp('Starting direct solver...');
tic
ds = DirectSolver( m, mu, ld, f1, f2, mfu, mfd, mfp, mim, ...
                   dirBoundID, dirBound1, dirBound2, ...
                   neuBoundID, neuBound1, neuBound2,...
                   'H1Semi');
Zt = ds.solve();
toc;

% load('example1_zBig.mat', 'zBig');

% zInterp = gf_compute(mfu0, zBig', 'interpolate on', mfu);
% %zExtrap = gf_compute(mfu, zBig', 'extrapolate on', mfu0);
% zInterp = zInterp';
% zInterp = ds.Q * (ds.Q' * zInterp) + ds.U0';
