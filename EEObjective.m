classdef EEObjective < handle
    properties
        F
        Z
        eps = 1E-6; 
        ds
        L
        K
        M
        KM
        BTpF
        BzCpTMIBzCp
        BzCp
    end

    methods
        function obj = EEObjective(ds, Z, eps)
            obj.F = ds.F;
            obj.Z = Z;

            f = ds.Q' * ds.F(1:ds.udof);
            zp = Z(1:ds.udof);

            z = ds.Q' * zp;
            p = Z( (ds.udof + 1):size(Z,1));

            obj.eps = eps;
            obj.L = ds.Q' * ds.getAdjointStiffness( zp );

            % K = 2*gf_asm( 'volumic',...
            % [ 't=comp(vGrad(#1).vGrad(#1));',...
            %   'e=(t{:,2,3,:,5,6}+t{:,3,2,:,5,6}',...
            %   '+ t{:,2,3,:,6,5}+t{:,3,2,:,6,5})/4;',...
            %   'M(#1,#1)+=sym(e(:,i,j,:,i,j))'], ds.mim, ds.mfu);
            
            K = gf_asm( 'volumic', ['t=comp(vGrad(#1).vGrad(#1));',...
                             'M(#1,#1)+=sym(t(:,i,j,:,i,j))'],ds.mim,ds.mfu);
                       
            % ds.K1Q = gf_asm('laplacian',ds.mim, ds.mfp, ds.mfp, gf_mesh_fem_get(ds.mfp, 'eval', {'1'} ));
            
            % ds.M1Q = gf_asm('mass matrix',ds.mim, ds.mfp, ds.mfp);

            %gf_compute(ds.mfu, z', 'H1 norm', ds.mim)

            
            % K1 =  2*gf_asm( 'volumic',...
            %    [ 'm=data$1(#2);',...
            %    't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
            %    'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
            %    ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/4;',...
            %    'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k).m(k));' ], ...
            %    ds.mim, ds.mfu, ds.mfd, ...
            %    gf_mesh_fem_get(ds.mfd, 'eval', {1} ) );

            % K2 = 2*gf_asm( 'volumic', ...
            %   [ 'm=data(#2);', ...
            %     't=comp(vGrad(#1).vGrad(#1).Base(#2));', ...
            %     's=(t(:,i,j,:,i,j,:) + t(:,i,j,:,j,i,:)', ...
            %     '+ t(:,j,i,:,i,j,:) + t(:,j,i,:,j,i,:))/4;', ...
            %     'M(#1,#1)+= sym(s(:,:,k).m(k))'], ...
            %     ds.mim, ds.mfu, ds.mfd,...
            %     gf_mesh_fem_get(ds.mfd, 'eval', {1} ) );
    

            % K3 = gf_asm( 'linear_elasticity', ds.mim, ds.mfu, ds.mfd, ...
            %             gf_mesh_fem_get(ds.mfd, 'eval', {0}),...
            %             gf_mesh_fem_get(ds.mfd, 'eval', {1}));

            % K4 = gf_asm('volumic', ...
            %       ['a=data$1(#2);', ...
            %        'M(#1,#1)+=sym(comp(vGrad(#1).vGrad(#1).Base(#2))(:,j,k,:,j,k,p).a(p))'],...
            %        ds.mim, ds.mfu, ds.mfd,...
            %        gf_mesh_fem_get(ds.mfd, 'eval', {1} ) );


            %keyboard
            M = ds.Mp;

            obj.KM = ds.Q' *(K + M ) * ds.Q;

            obj.ds = ds;
            Mp = gf_asm('mass matrix', ds.mim, ds.mfp);

            obj.BTpF = ds.B * p - f;
            obj.BzCp = (ds.B' * z - ds.C * p);
            obj.BzCpTMIBzCp = obj.BzCp' * (Mp \ obj.BzCp);


        end

        function [Fa, Ga] = evaluate(obj, A)

            LABTpF = obj.L * A  + obj.BTpF;

            KMILABTpF = obj.KM \ LABTpF;
            
            Fa = 0.5*LABTpF' * KMILABTpF; % + 0.5*obj.BzCpTMIBzCp;

            Fa = Fa + 0.5*obj.eps * A' * obj.ds.R * A;

            Ga = obj.L' *  KMILABTpF;


            Ga = Ga + obj.eps *  obj.ds.R * A;


            % Ha = obj.L' * obj.KMI * obj.L;

            % Ha = Ha + obj.eps *  obj.ds.K1 ;
        end 
    end
end
