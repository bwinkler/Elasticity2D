classdef OLSObjective < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        gradMeth
        M
        Mp
    end

    methods
        function obj = OLSObjective( ds, Z, eps, gradMeth)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            obj.gradMeth = gradMeth;
            if strcmp(gradMeth, 'Classical')
              %error('OLS classical gradient non-functional')
              obj.Kphi = ds.getKphi();
            end
            obj.M = ds.M;
            obj.Mp =  gf_asm('mass_matrix', ds.mim, ds.mfp);

        end

        function [Fa, Ga, Ha] = evaluate(obj, A)

            ds = obj.ds;
            ds.setParam( A );

            U = ds.solve();

            up = U(1:ds.udof);
            zp = obj.Z(1:ds.udof);
            
            % u = (ds.Q* (ds.Q' * up)) + ds.U0';
            % z = (ds.Q*(ds.Q' * zp)) + ds.U0';

            u = ds.Q' * up;
            z = ds.Q' * zp;

            pu = U( (ds.udof + 1): size(U,1));
            pz = obj.Z( (ds.udof + 1):size(U,1));


            Fa = 0.5 * ( u - z )' * ds.M * (u - z);

            %Fa = 0.5 *( (up-zp)' * ds.Mp * (up - zp));

            Fa = Fa + 0.5 * obj.eps * A' * ds.R * A;

            if strcmp(obj.gradMeth, 'Classical')
              Ga = zeros(ds.dof, 1);
               for i = [1:ds.dof]
                size(ds.Q'*obj.Kphi{i}*ds.Q )
                %size(u)
                du = A \  -(ds.Q'* obj.Kphi{i}*ds.Q  * u);
                Ga(i) = du'  * ds.M * (u-z);
              end
            elseif strcmp(obj.gradMeth, 'Hybrid')
              %keyboard
              L = ds.Q' * ds.getAdjointStiffness( up );
              k = size(ds.A,1);
              m = size(L,2);
              nn = size(ds.C,1);
              DUb = zeros(k,m);
              for i =1:m
                E = zeros(m,1);
                E(i) = 1;

                DU = ds.K \ [-L*E; zeros(nn,1)];
                DUb(:,i) = DU(1:k);
              end

              Ga = (u-z)' * ds.M *DUb;
              Ga = Ga';


              Ha = DUb'*ds.M * DUb;

              % W = ds.Kp \ [  zp - up; zeros(size(pu))]; 
              % Wb = W(1:ds.udof);
              % Lw = ds.getAdjointStiffness(Wb);
              % Ha = Ha + 2 * DUb' * ds.Q'*Lw;
              W = ds.K \ [ z - u; zeros(size(pu))];
              Wb = W(1:size(ds.A,1));
              Lw = ds.Q'* ds.getAdjointStiffness(ds.Q*Wb);
              Ha = Ha + 2 * DUb' * Lw;
              %norm( ds.A * Wb - Lw* A)
              Ha = 0.5 * (Ha + Ha');
              Ha = Ha + obj.eps * ds.R;

            elseif strcmp(obj.gradMeth, 'Adjoint')
              % W = ds.Kp \ [  zp - up; zeros(size(ds.C,1),1)]; 
              % Wb = W(1:size(ds.Ap,1));
              % L = ds.getAdjointStiffness(Wb);
              % Ga = up' * L;

              W = ds.K \ [ z - u; zeros(size(pu))];
              Wb = W(1:size(ds.A,1));
              L = ds.Q'* ds.getAdjointStiffness(ds.Q*Wb);
              Ga = u' * L ;

              %norm( ds.A * Wb - L* A)
              Ga = Ga';
            else
              L = ds.Q' * ds.getAdjointStiffness( up );
              Ga = -0.5*L'*(ds.A \  (ds.M * (u - z)));
              % L = ds.getAdjointStiffness( up );
              % Ga = -0.5 * L' * (ds.Ap \ (ds.Mp * (up-zp)));
              %Ga = -0.5*L'* inv(ds.A) * ds.M * (u - z);
            end
            Ga = Ga + obj.eps *  ds.R * A;
        end
    end
end
