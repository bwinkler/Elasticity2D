function dispMuComparison(ds, muExact, muComp, fileName, zLim)

  resolution = 300;
  format = 'png';
  errScale = 10;

  resString = sprintf('-r%d', resolution);
  fmtString = sprintf('-d%s', format);

  muDiff = muExact - muComp;

  figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  zlim( zLim );
  print ('-zbuffer', resString, fmtString, sprintf('%s_exact.%s', fileName, format));

  figure('Renderer', 'zbuffer')
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  zlim( zLim );
  print ('-zbuffer', resString, fmtString, sprintf('%s_comp.%s', fileName, format));

  figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, (muExact - muComp), 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  print ('-zbuffer', resString , fmtString, sprintf('%s_err.%s', fileName, format));

  figure;
  semilogy( abs( muDiff ) );
  xlim([ 0 length(muDiff) ] );
  print('-painters', resString, fmtString, sprintf('%s_errnodal.%s', fileName, format));

end
