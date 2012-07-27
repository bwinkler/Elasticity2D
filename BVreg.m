function [fx,gx] = BVreg(ds, A, e, edges)

            % BV Regularization
            
  %e = 0.00001;
  %E = gf_mesh_get(ds.mesh, 'edges', 'merge');

  fx = 0;
  gx = zeros(size(A));

  for i = [ 1:size(edges,2) ]
              pid1 = edges(1,i);
              pid2 = edges(2,i);
              x = A(pid1) - A(pid2);
              p1 = gf_mesh_fem_get(ds.mfd, 'basic dof nodes', pid1);
              p2 = gf_mesh_fem_get(ds.mfd, 'basic dof nodes', pid2);
              len = norm( p1 - p2 );
              f = sqrt( x^2 + e ) * len;
              %df =  x / f;
              df=(x.^3+2*e*x)./ ( x^2 + e ) .^(1.5);
              fx = fx + f;
              gx(pid1) = gx(pid1) + df * len;
              gx(pid2) = gx(pid2) - df * len;
            end 
end
