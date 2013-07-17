clear all;
gf_workspace('clear all');
warning off;

% Create the mesh
n = 73;
h = 1/(n+1);


m = gf_mesh('cartesian', [0:h:1], [0:h:1] );
%m = gf_mesh('triangles grid', [0:h:1], [0:h:1] );

pts = gf_mesh_get(m,'pts');
pidleft = find( abs(pts(1,:) ) < 1/(2*(n+1)) );
pidright = find( abs(pts(1,:) - 1 ) < 1/(2*(n+1)) );
pidbottom = find( abs(pts(2,:) ) < 1/(2*(n+1)) );
pidtop = find( abs(pts(2,:) - 1 ) < 1/(2*(n+1)) );

dirBoundFaces = gf_mesh_get(m, 'faces_from_pid', pidtop );
neuBoundFaces = gf_mesh_get(m, 'faces_from_pid',...
                  union( pidleft, union(pidright, pidbottom) ) );

dirBoundID = 1;
neuBoundID = 2;

gf_mesh_set(m, 'boundary', dirBoundID, dirBoundFaces);
gf_mesh_set(m, 'boundary', neuBoundID, neuBoundFaces);

mfu = gf_mesh_fem(m,2);
gf_mesh_fem_set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));

mfd = gf_mesh_fem(m);
gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_QK(2,1)'));

mfp = gf_mesh_fem(m);
gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_QK(2,0)'));

mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,2)'));
%mim = gf_mesh_im(m, gf_integ('IM_QUAD(2)'));

% mfu = gf_mesh_fem(m,2);
% gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_PK(2,2)'));

% mfd = gf_mesh_fem(m);
% gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_PK(2,1)'));

% mfp = gf_mesh_fem(m);
% gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_PK(2,0)'));
% mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(2)'));


%mu = '2.5+0.25*sin(2*pi*x)';
mu = '0.2*sin(pi*x).* ( x >= 0.2 & x <= 0.4 & y >=0.2 & y <= 0.4)';
mu = strcat(mu, '+ 0.5*sin(pi*x).* ( x >= 0.6 & x <= 0.8 & y >=0.6 & y <= 0.8)');
mu = strcat(mu, '+ 0.3*sin(pi*x).* ( x >= 0.2 & x <= 0.4 & y >= 0.6 & y <=0.8)');
mu = strcat(mu, '+ 0.1*sin(pi*x).* ( x >= 0.6 & x <= 0.8 & y >= 0.2 & y <=0.4)');
mu = strcat(mu, '+ 1');
ld = 1E6;

% For totally incompressible case.
%ld = 'Incompressible';

f1 = '1+0.1*x.*x';
f2 = '0.1*y';


dirBound1 = '0';
dirBound2 = '0';

neuBound1 = '0.5*1+x.^2';
neuBound2 = '0';


disp('Starting direct solver...');
tic
ds = DirectSolver( m, mu, ld, f1, f2, mfu, mfd, mfp, mim, ...
                   dirBoundID, dirBound1, dirBound2, ...
                   neuBoundID, neuBound1, neuBound2,...
                   'BV');
Zt = ds.solve();
toc;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

pEst = estimateP( ds, zExact ); 

%Z = [ zExact; zeros(size(pExact))];
Z = [ zExact; pEst ];

%Z = Zt;
