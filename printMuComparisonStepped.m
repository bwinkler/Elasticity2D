function printMuComparisonStepped(varargin)
%function printMuComparisonStepped(ds, muExact, muComp, muHist, fileName, zLim, stepZLim)

  ds = varargin{1};
  muExact = varargin{2};
  muComp  = varargin{3};
  muHist = varargin{4};
  fileName = varargin{5};
  zLim = varargin{6};
  stepZLim = varargin{7};

  if(nargin > 7)
    optSteps = varargin{8};
  else
    optSteps = 1:size(muHist,1);
  end

  if(nargin > 8)
    rotation = varargin{9};
  else
    rotation = [-37.5, 30 ];
  end

  resolution = 600;
  format = 'eps';
  errScale = 0.0000000000005;

  resString = sprintf('-r%d', resolution);
  fmtString = sprintf('-d%s', format);

  muDiff = muExact - muComp;

  %style = hgexport('readstyle','SIAM02');

  figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  zlim( zLim );
  view(rotation);
  %print ('-zbuffer', resString, fmtString, sprintf('%s_exact.%s', fileName, format));
  %hgexport(gcf,sprintf('%s_exact.%s', fileName, format), style);
  savefig(sprintf('%s_exact', fileName), 'eps', '-rgb', '-r600', '-crop');

  figure('Renderer', 'zbuffer')
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  zlim( zLim );
  view(rotation);
  %print ('-zbuffer', resString, fmtString, sprintf('%s_comp.%s', fileName, format));
  %hgexport(gcf,sprintf('%s_comp.%s', fileName, format), style);
  savefig(sprintf('%s_comp', fileName), 'eps', '-rgb', '-r600', '-crop');

  for i = optSteps
    figure('Renderer', 'zbuffer')
    gf_plot( ds.mfd, muHist(i,:),  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
    zlim( stepZLim );
    view(rotation);
    %hgexport(gcf,sprintf('%s_comp_%d.%s', fileName, i, format), style);
    %print ('-zbuffer', resString, fmtString, sprintf('%s_comp_%d.%s', fileName, i, format));
    savefig(sprintf('%s_comp_%d', fileName,i), 'eps', '-rgb', '-r600', '-crop');
  end

  figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, abs(muDiff), 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off', 'deformation_scale', 1 );
  set(gca, 'ZScale', 'log');
  view(rotation);
  %hgexport(gcf,sprintf('%s_err.%s', fileName, format), style);
  %print ('-zbuffer', resString , fmtString, sprintf('%s_err.%s', fileName, format));
  savefig(sprintf('%s_err', fileName), 'eps', '-rgb', '-r600', '-crop');

  % figure;
  % semilogy( abs( muDiff ) );
  % xlim([ 0 length(muDiff) ] );
  % print('-painters', resString, fmtString, sprintf('%s_errnodal.%s', fileName, format));

end
