classdef MOLSObjectiveJ1 < handle
    properties
        ds
        Z
        Kphi
        eps = 1E-6; 
        gradMeth
        edges
        BtCiB
    end

    methods
        function obj = MOLSObjectiveJ1( ds, Z, eps, gradMeth)
            obj.ds = ds;
            obj.Z = Z;
            obj.eps = eps;
            obj.gradMeth = gradMeth;

            % These both should likely be move to DirectSolver
            if strcmp(obj.gradMeth, 'Classical')
              obj.Kphi = ds.getKphi();
            end
            % if strcmp(obj.gradMeth, 'Adjoint')
            %   obj.Kphi = ds.getKphi();
            % end
            if strcmp(ds.regType, 'BV')
              obj.edges = gf_mesh_get(ds.mesh, 'edges', 'merge'); 
            end

            obj.BtCiB =  ds.Bp* ( ds.C \ ds.Bp');
            %obj.BtCiB =  (ds.Bp / ds.C) * ds.Bp';

            %obj.BtCiB = obj.BtCiB +  1E-4*speye(size(obj.BtCiB));
            %obj.BtCiB = ds.Bp * inv(ds.C) * ds.Bp';
            %obj.BtCiB = ds.B *  (ds.C \ ds.B');
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
            Fa = 0.5 * ( (up - zp)' * ds.Ap * (up - zp)... 
                         + (pu - pz)' * ds.C * (pu-pz) );

            % Fa = 0.5 * ( (u - z)' * ds.A * (u - z)... 
            %              + (pu - pz)' * ds.C * (pu-pz) );
            % Second trilinear form
            % Fa = 0.5 * ( (up - zp)' * ds.Ap * (up - zp)... 
            %               - (pu - pz)' * ds.C * (pu-pz) ) ...
            %               + (up -zp)' * ds.Bp * (pu-pz);


            if strcmp(obj.gradMeth, 'Classical')
              % Classical gradient
              Ga = zeros(ds.dof, 1);
              for i = [1:ds.dof]
                Ga(i) = -0.5 * (up + zp)' * obj.Kphi{i} * (up - zp);
              end
            elseif strcmp(obj.gradMeth, 'Adjoint')
              %keyboard
              % W = ds.K \ [-ds.A*(u-z) - ds.B*(pu-pz); zeros(size(pu))];
              % wp = ds.Q*W(1:size(u,1));
              % Luz = ds.Q'*ds.getAdjointStiffness(up-zp);
              % Lw = ds.Q'*ds.getAdjointStiffness(wp);
              % Ga =0.5 * (u-z)' * Luz + u'*Lw;
              
              W = ds.Kp \ [-ds.Ap*(up-zp) - ds.Bp*(pu-pz); zeros(size(pu))];
              wp = W(1:size(up,1));
              Luz = ds.getAdjointStiffness(up-zp);
              Lw = ds.getAdjointStiffness(wp);
              Ga = 0.5* (up-zp)' *Luz + up' * Lw;
              Ga = Ga';
            else
              % Adjoint-Stiffness Gradient 
              
              % J1

              subMeth = 2;

              if subMeth == 1
                Lu = ds.getAdjointStiffness(up);
                Luz = ds.getAdjointStiffness(up-zp);
                BtCiBLui = (ds.Ap + obj.BtCiB) \ Lu;
                Ga = -0.5*(up + zp)' * Luz;
                Ga = Ga - (pu - pz)' * ds.Bp' * BtCiBLui;
                Ga = Ga + (up - zp)' * obj.BtCiB * BtCiBLui;
              elseif subMeth == 2
                Lu = ds.getAdjointStiffness(up);
                Luz = ds.getAdjointStiffness(up-zp);
                BtCiBLui = (ds.Ap + obj.BtCiB) \ Lu;
                Ga = -0.5*(up + zp)' * Luz;
                Ga = Ga - (pu - pz)' * ds.Bp' * BtCiBLui;
                BpTApI = ds.Bp' / ds.Ap;
                Ga = Ga + (up-zp)' * ds.Bp * ( (ds.C + BpTApI * ds.Bp ) \ (BpTApI * Lu));
                %Ga = Ga + (up-zp)' * ds.Bp * ( (-ds.C - ds.Bp'* (ds.Ap \ ds.Bp)) \ ( ds.Bp' * (ds.Ap \ Lu )));
              elseif subMeth == 3
                Lu = ds.getAdjointStiffness(up);
                Luz = ds.getAdjointStiffness(up-zp);
                Ga = -0.5*(up + zp)' * Luz;

                % n = size(up,1);
                % k = size(pu,1);
                % V = zeros(size(ds.Kp,1), n);
                % for i = 1:n
                %   keyboard
                %   e = spalloc(size(ds.Kp,1), 1,1);
                %   e(i) = 1;
                %   V(:,i) = ds.Kp \ e;
                % end

                n = size(u,1);
                k = size(pu,1);
                V = zeros(n+k, n);
                for i = 1:n
                  e = spalloc(n+k, 1,1);
                  e(i) = 1;
                  V(:,i) = ds.K \ e;
                end
                % Ga = Ga - (pu - pz)' * ds.Bp' * V(1:n,1:n) * Lu;
                % Ga = Ga + (up-zp)' * ds.Bp * V(n+1:n+k,1:n) * Lu;

                Ga = Ga - (pu - pz)' * ds.B' * V(1:n,1:n) * ds.Q'*Lu;
                Ga = Ga + (u-z)' * ds.B * V(n+1:n+k,1:n) * ds.Q'*Lu;
              end
              Ga = Ga';


            end

            if strcmp(ds.regType, 'BV')
              % Use BV regularization.
              [fx, gx] = BVreg(ds, A, 0.00001, obj.edges);
              Fa = Fa + obj.eps * fx;
              Ga = Ga + obj.eps * gx;
            else
              % Use smooth regularization.
              Ga = Ga + obj.eps *  ds.R * A;
              Fa = Fa +  0.5*obj.eps * A' * ds.R * A;
            end
            
        end
    end
end
