function dispMuComparison( varargin )
%function dispMuComparison(ds, muExact, muComp)

  ds = varargin{1};
  muExact = varargin{2};
  muComp = varargin{3};
  if(nargin == 4)
    rotation = varargin{4};
  else
    rotation = [-37.5, 30 ];
  end
  figure('Renderer', 'zbuffer');
  subplot(1,2,1);
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  title('Exact \mu');
  view(rotation);
  subplot(1,2,2);
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  title('Computed \mu');
  view(rotation);

  disp( sprintf('Error 2-norm: %f', norm( muExact - muComp ) ));
  disp( sprintf('Error Inf-norm: %f', norm(muExact - muComp, Inf) ) );

  M = gf_asm('mass_matrix', ds.mim, ds.mfd);
  K = gf_asm('laplacian', ds.mim, ds.mfd, ds.mfd,gf_mesh_fem_get(ds.mfd, 'eval', {1} ) );

  disp( sprintf('Error L2-norm: %g', (muExact-muComp)* M * (muExact-muComp)'));
  disp( sprintf('Error H1 semi-norm: %g', (muExact-muComp)* K * (muExact-muComp)'));
  disp( sprintf('Error H1 norm: %g', (muExact-muComp)* (K+M) * (muExact-muComp)'));
end
