clear all;
gf_workspace('clear all');
warning off;

% Create the mesh
n = 20;
h = 1/(n+1);
m = gf_mesh('cartesian', [0:h:1], [0:h:1] );


pts = gf_mesh_get(m,'pts');
pidleft = find( abs(pts(1,:) ) < 1/(2*(n+1)) );
pidright = find( abs(pts(1,:) - 1 ) < 1/(2*(n+1)) );
pidbottom = find( abs(pts(2,:) ) < 1/(2*(n+1)) );
pidtop = find( abs(pts(2,:) - 1 ) < 1/(2*(n+1)) );

dirBoundFaces = gf_mesh_get(m, 'faces_from_pid', union( pidleft, pidbottom) );
neuBoundFaces = gf_mesh_get(m, 'faces_from_pid', union( pidright, pidtop) );

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

mu = '1+x.*y.*sin(pi*x).*sin(pi*y)';

ld = 1E6;

f2 = 'cos(pi*x)';
f1 = '-0.2*x';

% dirBound1 = '-0.01*sin(pi*y)';
% dirBound2 = '0.1*sin(x )';

dirBound1 = '0.1*sin(pi*y)';
dirBound2 = '0.1*sin(pi*x)';

% dirBound1 = '0';
% dirBound2 = '0';

neuBound1 = '0.1+x';
neuBound2 = '0.1+y';

tic
ds = DirectSolver( m, mu, ld, f1, f2, mfu, mfd, mfp, mim,... 
                   dirBoundID, dirBound1, dirBound2, ...
                   neuBoundID, neuBound1, neuBound2, ...
                   'H1Semi');
Zt = ds.solve();

toc;

zExact = Zt(1:ds.udof);
pExact = Zt(ds.udof + 1: length(Zt));


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



