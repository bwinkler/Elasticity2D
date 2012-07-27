clear all;
gf_workspace('clear all');
warning off;

% Create the mesh
n = 50;
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

mfu = gf_mesh_fem(m,2);
gf_mesh_fem_set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));

mfd = gf_mesh_fem(m);
gf_mesh_fem_set(mfd, 'fem', gf_fem('FEM_QK(2,1)'));

mfp = gf_mesh_fem(m);
gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_QK(2,0)'));

mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,2)'));

mu = '1+x.*y.*sin(pi*x).*sin(pi*y)';

ld = 1E6;

%mu = '2+x.^2+y.^2';
fx = '2+sin(2*pi*x)';
fy = '2+0.1*y';

dirBoundx = 'sin(x)';
dirBoundy = '0';

neuBoundx = '1';
neuBoundy = '1';

tic
ds = DirectSolver( m, mu, ld, fx, fy, mfu, mfd, mfp, mim,... 
                   TOP, dirBoundx, dirBoundy, ...
                   REST, neuBoundx, neuBoundy, ...
                   'H1Semi');
Zt = ds.solve();

toc;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));

Z = [ zExact; zeros(size(pExact))];


% u1 = Ut(1:ds.udof);
% p = Ut((ds.udof+1):length(Ut));


% udx = u1(1:2:end);
% udy = u1(2:2:end);

% noisep = 0.01;

% [udx, udy, err_norm] = addNoise( udx, udy, noisep, 'uniform');

% ud1 = [ udx, udy ]';

% M = gf_asm('mass matrix',mim, mfu);
% alpha = 1E-9;
% K = alpha * gf_asm( 'volumic',...
%              [... 
%              't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
%              'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
%              ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
%              'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k));' ], ...
%              ds.mim, ds.mfu, ds.mfd );
% F = gf_asm( 'volumic_source', mim, mfu, mfd, ud1);


% %QB = gf_asm('boundary qu term', BORDER, mim, mfu, mfd, [zeros(1,ds.dof); zeros(1,ds.dof); zeros(1,ds.dof); zeros(1,ds.dof);  ]);

% z1 = ds.Q * ((ds.Q' * K * ds.Q + ds.Q' * M * ds.Q) \ (ds.Q'*F)) + ds.U0';

% zd = [];


% for i = [1:size(udx)]
%    zd = [ zd; udx(i); udy(i) ];
% end

% zd = ds.Q * (ds.Q' * zd) + ds.U0';

% deform_scale = 0.1;
% figure;

% subplot(2,1,1);
% gf_plot( ds.mfu, u1', ...
%     'deformed_mesh', 'on',...
%     'norm', 'on', ...
%     'deformation', u1', ...
%     'deformation_mf', ds.mfu, ...
%     'deformation_scale', deform_scale, ...
%     'disp_options', 'off', ... 
%     'title', 'Displacement' );
% title('Displacement and Deformation');

% subplot(2,1,2);
% gf_plot( ds.mfu, z1', ...
%     'norm', 'on', ...
%     'deformed_mesh', 'on',...
%     'quiver', 'off', ...
%     'deformation', z1', ...
%     'deformation_mf', ds.mfu, ...
%     'deformation_scale', deform_scale, ...
%     'disp_options', 'off', ... 
%     'title', 'Displacement' );
% title('Displacement and Deformation smoothed');



% % %colorbar;

%  return



