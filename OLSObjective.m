classdef OLSObjective < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        args
    end

    methods
        function obj = OLSObjective( ds, Z, eps)
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
            
            L = ds.Q' * ds.getAdjointStiffness( u );
            u = ds.Q' * u;
            z = ds.Q' * obj.Z(1:ds.udof);

            %Fa = 0.5 * (U - obj.Z)' * ds.M * (U - obj.Z);

            Fa = 0.5 * ( u - z )' * ds.M * (u - z);
            Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;

            %Ga = zeros(ds.dof, 1);


            Ga = -0.5 * L' * inv( ds.A ) * ds.M * (u - z);

            %for i = [1:ds.dof]
            %  %dU = ds.Kp \ (-obj.Kphi{i} * U);
            %  du = ds.Ap \ (-obj.Kphi{i}(1:ds.udof,1:ds.udof) * u);
            %  %dU =gf_linsolve('superlu', ds.Kp, -obj.Kphi{i} * U );
            %  
            %  %du =gf_linsolve('gmres', ds.Ap, (-obj.Kphi{i}(1:ds.udof,1:ds.udof) * u));

            %  %L =  -( ds.Ap \ (ds.M * (u - z) ));
            %  %Ga(i) = (U - obj.Z)' * ds.M * dU;
            %  Ga(i) = ( u - z )' * ds.M * du;
            %  %
            %  %Ga(i) = u' * obj.Kphi{i}(1:ds.udof,1:ds.udof)  * L;
            %end

            Ga = Ga + 0.5 * obj.eps *  ds.R * A;
        end
    end
end
