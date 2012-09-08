function uSmooth = smoothData( u, alpha, ds )

  [u1, u2] = splitVector(u);

  Ms = gf_asm('mass matrix', ds.mim, ds.mfd);
  Ks = alpha * gf_asm('laplacian', ds.mim, ds.mfd, ds.mfd,...
                       gf_mesh_fem_get(ds.mfd, 'eval', {1}));

  z1 = gf_asm('volumic source', ds.mim, ds.mfd,ds.mfd, u1');
  z2 = gf_asm('volumic source', ds.mim, ds.mfd, ds.mfd, u2');

  u1s = (Ks + Ms) \ z1;
  u2s = (Ks + Ms) \ z2;

  uSmooth = mergeVectors(u1, u2s);

  uSmooth = ds.Q * (ds.Q' * uSmooth) + ds.U0';

end
