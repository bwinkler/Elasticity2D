function dispDisplacementComparison(ds, u, z, deformScale)
  figure('Renderer', 'zbuffer');

  subplot(2,1,1);
  gf_plot( ds.mfu, u, ...
           'deformed_mesh', 'on',...
           'norm', 'on', ...
           'deformation', u, ...
           'deformation_mf', ds.mfu, ...
           'deformation_scale', deformScale, ...
           'disp_options', 'off' );
  title('Displacement u');

  subplot(2,1,2);
  gf_plot( ds.mfu, z, ...
           'norm', 'on', ...
           'deformed_mesh', 'on',...
           'deformation', z, ...
           'deformation_mf', ds.mfu, ...
           'deformation_scale', deformScale, ...
           'disp_options', 'off');
  title('Deformation z');

end
