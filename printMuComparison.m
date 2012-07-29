function dispMuComparison(ds, muExact, muComp, fileName, zLim)
  figure('Renderer', 'zbuffer');
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  zlim( zLim );
  print ('-zbuffer', '-r600', '-dtiff', sprintf('%s_exact.tiff', fileName));
  figure('Renderer', 'zbuffer')
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  zlim( zLim );
  print ('-zbuffer', '-r600', '-dtiff', sprintf('%s_comp.tiff', fileName));

end
