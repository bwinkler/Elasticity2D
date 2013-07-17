classdef MOLSObjectiveJ2 < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        gradMeth
        edges
    end

    methods
        function obj = MOLSObjectiveJ2( ds, Z, eps, gradMeth)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            obj.gradMeth = gradMeth;

            % These both should likely be move to DirectSolver
            if strcmp(gradMeth, 'Classical')
              obj.Kphi = ds.getKphi();
            end
            if strcmp(ds.regType, 'BV')
              obj.edges = gf_mesh_get(ds.mesh, 'edges', 'merge'); 
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
            
            pu = U( (ds.udof + 1): size(U,1));
            pz = obj.Z( (ds.udof + 1):size(U,1));


            % First trilinear form
            % Fa = 0.5 * ( (up - zp)' * ds.Ap * (up - zp)... 
            %              + (pu - pz)' * ds.C * (pu-pz) );

            % Second trilinear form
            % Fa = 0.5 * ( (up - zp)' * ds.Ap * (up - zp)... 
            %               - (pu - pz)' * ds.C * (pu-pz) ) ...
            %               + (up -zp)' * ds.Bp * (pu-pz);

            Fa = 0.5 * ( (u - z)' * ds.A * (u - z)... 
                           - (pu - pz)' * ds.C * (pu-pz) ) ...
                           + (u -z)' * ds.B * (pu-pz);

            if strcmp(obj.gradMeth, 'Classical')
              % Classical gradient
              Ga = zeros(ds.dof, 1);

              for i = [1:ds.dof]
                Ga(i) = -0.5 * (up + zp)' * obj.Kphi{i} * (up - zp);
              end
            else
              % Adjoint-Stiffness Gradient 
              
              % Second trilinear form
              % L = ds.getAdjointStiffness(up + zp); 
              % Ga = -0.5* L' * (up - zp);
              Lu = ds.getAdjointStiffness(up);
              Lz = ds.getAdjointStiffness(zp);
              Ga = -0.5 * Lu' * up + 0.5 * Lz' * zp;
            end

            if strcmp(ds.regType, 'BV')
              % Use BV regularization.
              [fx, gx] = BVreg(ds, A, 0.00001, obj.edges);
              Fa = Fa + obj.eps * fx;
              Ga = Ga + obj.eps * gx;
            else
              % Use smooth regularization.
              Ga = Ga + obj.eps *  ds.R * A;
              Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;
            end
            
        end
    end
end
