classdef OLSObjective < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        gradMeth
    end

    methods
        function obj = OLSObjective( ds, Z, eps, gradMeth)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            obj.gradMeth = gradMeth;
            if strcmp(gradMeth, 'Classical')
              error('OLS classical gradient non-functional')
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
            pu = U( (ds.udof + 1): size(U,1));
            pz = obj.Z( (ds.udof + 1):size(U,1));


            %Fa = 0.5 * ( u - z )' * ds.M * (u - z);
            Fa = 0.5 * ( (up - zp)' * ds.Mp * (up - zp)... 
                          - (pu - pz)' * ds.C * (pu-pz) ) ...
                          + (up -zp)' * ds.Bp * (pu-pz);

            Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;

            if strcmp(obj.gradMeth, 'Classical')
              Ga = zeros(ds.dof, 1);
               for i = [1:ds.dof]
                du = ds.A /  -(( ds.Q'*obj.Kphi{i}*ds.Q)  * u);
                Ga(i) = du'  * ds.M * (u-z);
              end
            else
              L = ds.Q' * ds.getAdjointStiffness( up );
              Ga = -0.5*L' *  (ds.A \  (ds.M * (u - z)));
            end
            Ga = Ga + 0.5* obj.eps *  ds.R * A;
        end
    end
end
