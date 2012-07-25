function dispMuComparison(ds, muExact, muComp)
  figure('Renderer', 'zbuffer');
  subplot(1,2,1);
  gf_plot( ds.mfd, muExact, 'zplot', 'on', 'deformed_mesh', 'on', 'disp_options','off');
  title('Exact \mu');
  subplot(1,2,2);
  gf_plot( ds.mfd, muComp,  'zplot', 'on', 'deformed_mesh', 'on', 'disp_options', 'off');
  title('Computed \mu');
end
