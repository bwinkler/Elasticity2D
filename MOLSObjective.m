classdef MOLSObjective < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        gradMeth
    end

    methods
        function obj = MOLSObjective( ds, Z, eps, gradMeth)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            obj.gradMeth = gradMeth;

            if strcmp(gradMeth, 'Classical')
              obj.Kphi = ds.getKphi();
            end

        end

        function [Fa, Ga] = evaluate(obj, A)

            ds = obj.ds;
            ds.setParam( A );

            U = ds.solve();
            up = U(1:ds.udof);
            zp = obj.Z(1:ds.udof);


            u = ds.Q' * up;
            z = ds.Q' * zp;


            % BV Regularization
            
            % e = 0.00001;
            % E = gf_mesh_get(ds.mesh, 'edges', 'merge');

            % fx = 0;
            % gx = zeros(size(A));

            % for i = [ 1:size(E,2) ]
            %   pid1 = E(1,i);
            %   pid2 = E(2,i);
            %   x = A(pid1) - A(pid2);
            %   p1 = gf_mesh_fem_get(ds.mfd, 'basic dof nodes', pid1);
            %   p2 = gf_mesh_fem_get(ds.mfd, 'basic dof nodes', pid2);
            %   len = norm( p1 - p2 );
            %   f = sqrt( x^2 + e ) * len;
            %   df =  x / f;
            %   %df=(x.^3+2*e*x)./ ( x^2 + e ) .^(1.5);
            %   fx = fx + f;
            %   gx(pid1) = gx(pid1) + df * len;
            %   gx(pid2) = gx(pid2) - df * len;
            % end 


            Fa = 0.5 * (u - z)' * ds.A * (u - z);

            % Tikhonov Regularization
            Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;
            
            % BV
            %Fa = Fa + obj.eps * fx;

            if strcmp(obj.gradMeth, 'Classical')
              % Classical gradient
              Ga = zeros(ds.dof, 1);

              for i = [1:ds.dof]
                Ga(i) = -0.5 * (up + zp)' * obj.Kphi{i} * (up - zp);
              end
            else
              % Adjoint-Stiffness Gradient
              L = ds.Q' * ds.getAdjointStiffness(up + zp); 
              Ga = -0.5* L' * (u - z);
            end

            % Tikhonov Regularization
            Ga = Ga +  obj.eps *  ds.R * A;

            % BV
            %Ga = Ga + obj.eps * gx;
            
        end
    end
end
