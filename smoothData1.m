function uSmooth = smoothData( u, alpha, ds )

  Mp = gf_asm('mass matrix',ds.mim, ds.mfu);

  Kp = alpha * gf_asm( 'volumic',...
             [... 
             't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
             'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
             ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
             'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k));' ], ...
             ds.mim, ds.mfu, ds.mfd );
  M = ds.Q' * Mp * ds.Q;
  K = ds.Q' * Kp * ds.Q;



  [u1, u2 ] = splitVector(u);

  ud = [u1,u2]';
  Z = gf_asm( 'volumic_source', ds.mim, ds.mfu, ds.mfd, ud);

  %uSmooth = ds.Q * ((ds.Q' * K * ds.Q + ds.Q' * M * ds.Q) \ (ds.Q'*F));
  uSmooth = (K+M) \ ( ds.Q' * ( Z - (Kp + Mp) * ds.U0') );
  uSmooth = ds.Q * uSmooth + ds.U0';


  % [u1, u2] = splitVector(u);

  % Ms = gf_asm('mass matrix', ds.mim, ds.mfd);
  % Ks = alpha * gf_asm('laplacian', ds.mim, ds.mfd, ds.mfd,...
  %                      gf_mesh_fem_get(ds.mfd, 'eval', {1}));

  % z1 = gf_asm('volumic source', ds.mim, ds.mfd,ds.mfd, u1');
  % z2 = gf_asm('volumic source', ds.mim, ds.mfd, ds.mfd, u2');

  % u1s = (Ks + Ms) \ z1;
  % u2s = (Ks + Ms) \ z2;

  % uSmooth = mergeVectors(u1s, u2s);

  % uSmooth = ds.Q * (ds.Q' * uSmooth) + ds.U0';




end
