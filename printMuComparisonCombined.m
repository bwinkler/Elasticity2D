function dispMuComparisonCombined(ds, muExact, muComp, fileName, zLim)

  resolution = 300;
  format = 'png';
  errScale = 10;

  resString = sprintf('-r%d', resolution);
  fmtString = sprintf('-d%s', format);

  muDiff = muExact - muComp;

  figure('Renderer', 'zbuffer');
  subplot(2,2,1);
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  zlim( zLim );
  title('Exact \mu');
  %print ('-zbuffer', resString, fmtString, sprintf('%s_exact.%s', fileName, format));

  %figure('Renderer', 'zbuffer')
  subplot(2,2,2);
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  zlim( zLim );
  title('Computed \mu');

  %print ('-zbuffer', resString, fmtString, sprintf('%s_comp.%s', fileName, format));

  subplot(2,2,3);
  %figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, (muExact - muComp), 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  title('Coefficient Error');
  %print ('-zbuffer', resString , fmtString, sprintf('%s_err.%s', fileName, format));

  %figure;
  subplot(2,2,4);
  semilogy( abs( muDiff ) );
  xlim([ 0 length(muDiff) ] );
  title('Nodal Error');
  print('-painters', resString, fmtString, sprintf('%s_combined.%s', fileName, format));

end
