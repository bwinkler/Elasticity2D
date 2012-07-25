classdef MOLSObjective < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        args
    end

    methods
        function obj = MOLSObjective( ds, Z, eps)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            %obj.Kphi = ds.getKphi();
        end

        function [Fa, Ga] = evaluate(obj, A)

            ds = obj.ds;
            ds.setParam( A );

            U = ds.solve();
            u = U(1:ds.udof);
            z = obj.Z(1:ds.udof);

            L = ds.Q' * ds.getAdjointStiffness(u + z); 

            u = ds.Q' * u;
            z = ds.Q' * z;


            %norm( L*A - ds.Ap * (u + z) )

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
            %Fa = 0.5 * (U - obj.Z)' * ds.Kp * (U - obj.Z);

            % Tikhonov Regularization
            Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;
            
            % BV
            %Fa = Fa + obj.eps * fx;

            %size(L')

            % Adjoint-Stiffness Gradient
            Ga = -0.5* L' * (u - z);

            % Classical gradient
            %Ga = zeros(ds.dof, 1);

            %for i = [1:ds.dof]
              %Ga(i) = -0.5 * (u + z)' * obj.Kphi{i} * (u - z);
            %  Ga(i) = -0.5 * (U + obj.Z)' * obj.Kphi{i} * (U - obj.Z);
            %end

            % Tikhonov Regularization
            Ga = Ga +  obj.eps *  ds.R * A;

            % BV
            %Ga = Ga + obj.eps * gx;
            
        end
    end
end
