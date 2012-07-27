classdef DirectSolver < handle
    properties
        n = 20;    % number of internal nodes
        mu         % The mu coefficient function
        fx         % load on x axis
        fy         % load on y axis
        mesh       % The mesh to solve on
        mfu        % The solution mesh fem 
        mfd        % The parameter mesh fem
        mfp        % The pressure mesh fem
        mim        % getfem integration method
        dof        % The number of degrees of freedom
        fVec       % The discretized f(x)
        mVec       % The discretized mu(x)
        dirVec     % Discretized Dirichlet conditions
        neuVec     % Discretized Neumann conditions
        dirBoundID    % The boundary ID for the Dirichlet conditions
        K          % The stiffness matrix
        F          % The discretized F after dirichlet conditions
        M          % Mass matrix
        Mp
        K1         % The (n+1)x(n+1) stiffness matrix
        A
        B
        Q          % Projector into Dirichlet solution space
        U0
        Kp          % List of matrices
        Bp
        udof
        pdof
        R          % Regularization
        C
        Ap
        Fp
        Fn
        regType    % Regularization type        
    end
    methods
        function ds = DirectSolver( mesh, mu, ld, fx, fy, mfu, mfd, mfp, mim,... 
                                    dirBoundID, dirBound_x, dirBound_y,...
                                    neuBoundID, neuBound_x, neuBound_y, ...
                                    regType )

            % Solves the direct incompressible problem on a 2D cartesian mesh
            %
            % Parameters:
            % n - the number of internal nodes
            % mu - Coefficient function as string or as n+2 vector of values
            % f - Load function as string or as n+2 vector of values

            % To be useful, this will need to also take the boundaries

            % Initialize solver valueb
            ds.mesh =mesh;
            ds.mu = mu;
            ds.fx = fx;
            ds.fx = fy;

            ds.mfu = mfu;
            ds.mfd = mfd;
            ds.mfp = mfp;

            ds.mim = mim;

            ds.regType = regType;

            ds.mVec = gf_mesh_fem_get(mfd, 'eval', {mu} );
            ds.fVec = gf_mesh_fem_get(mfd, 'eval', {fx;fy});

            ds.dirBoundID = dirBoundID;

            ds.dirVec = gf_mesh_fem_get(mfd, 'eval', { dirBound_x; dirBound_y } );

            ds.neuVec = gf_mesh_fem_get(mfd, 'eval', { neuBound_x; neuBound_y } );

            ds.Ap = gf_asm( 'volumic',...
                           [ 'm=data$1(#2);',...
                           't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
                           'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
                           ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
                           'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k).m(k));' ], ...
                           ds.mim, ds.mfu, ds.mfd, ds.mVec );
            ds.Bp = gf_asm( 'volumic', ...
                            'M(#1,#2)+=comp(vGrad(#1).Base(#2))(:,i,i,:)', ...
                            mim, mfu, mfp );

            ds.Fp = gf_asm( 'volumic_source', mim, mfu, mfd, ds.fVec );


            ds.Fn = gf_asm( 'boundary', neuBoundID, ...
                ['F=data(qdim(#1),#2);', ...
                    'V(#1)+=comp(vBase(#1).Base(#2))(:,i,j).F(i,j);'], ...
                    mim, mfu, mfd, ds.neuVec );

            
            if strcmp(ld, 'Incompressible')
              ds.C = zeros(size(ds.Bp,2));
            else

              ldInv = ld^(-1);
              ds.C = -gf_asm( 'volumic', ...
                          ['ld=data$1(#2);', ...
                            'M(#1,#1)+=comp(vBase(#1).vBase(#1).Base(#2))(:,i,:,i,j).ld(j)'], ...
                            mim, mfp, mfd, gf_mesh_fem_get(mfd, 'eval', ldInv ));
            end


            ds.Kp = [ds.Ap, ds.Bp; -ds.Bp', ds.C];
            
            ds.dof = gf_mesh_fem_get(mfd,  'nbdof');
            ds.udof = gf_mesh_fem_get(mfu, 'nbdof');
            ds.pdof = gf_mesh_fem_get(mfp, 'nbdof');

            [H, R] = gf_asm('dirichlet', dirBoundID, mim, mfu, mfd, ...
                            repmat( eye(2), [1,1,ds.dof] ) , ds.dirVec );
            [ds.Q, ds.U0] = gf_spmat_get(H, 'dirichlet_nullspace', R);

            ds.A = ds.Q' * ds.Ap * ds.Q;
            %Ft = ds.Q' * ds.Fp; 
            Ft = ds.Q' * ( (ds.Fp + ds.Fn) - ds.Ap * ds.U0');
            ds.B = ds.Q' * ds.Bp;


            G = zeros(size(ds.B,2), 1);

            ds.K = [ds.A, ds.B; ds.B', ds.C];
            ds.F = [Ft; G];

            ds.Mp = gf_asm('mass_matrix', mim, mfu);
            ds.M = ds.Q' * ds.Mp * ds.Q;


            switch ds.regType
            case 'L2'
              ds.R = gf_asm( 'mass_matrix', ds.mim, ds.mfd );
            case 'H1'
              ds.R = gf_asm( 'laplacian', ds.mim, ds.mfd, ds.mfd,...
                    gf_mesh_fem_get(mfd, 'eval', {1} ));
            case 'H1Semi'
              ds.R = gf_asm( 'mass_matrix', ds.mim, ds.mfd )...
                      + gf_asm( 'laplacian', ds.mim, ds.mfd, ds.mfd,...
                      gf_mesh_fem_get(mfd, 'eval', {1} ));
            case 'BV'
              ...
            otherwise
              error('Invalid regularization type.');
            end

          end
        function Ut = solve(ds)
          UP = ds.K \ ds.F;

          %UP =gf_linsolve('superlu', ds.K, ds.F );

          Ut = UP(1:size(ds.A,1), 1);
          P = UP((size(ds.A,1) + 1):size(UP,1));
          U = ds.Q*Ut + ds.U0';
          Ut = [U;P];

        end
        function setParam( ds, A )
          ds.mVec = A;

          ds.Ap = gf_asm( 'volumic',...
                           [ 'm=data$1(#2);',...
                           't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
                           'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
                           ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
                           'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k).m(k));' ], ...
                           ds.mim, ds.mfu, ds.mfd, ds.mVec );
          [H, R] = gf_asm('dirichlet', ds.dirBoundID, ds.mim, ds.mfu, ds.mfd, ...
                            repmat( eye(2), [1,1,ds.dof] ) , ds.dirVec );
          [ds.Q, ds.U0] = gf_spmat_get(H, 'dirichlet_nullspace', R);



          ds.A = ds.Q' * ds.Ap * ds.Q;
          ds.K = [ ds.A, ds.B; ds.B', ds.C ];
          ds.Kp = [ds.Ap, ds.Bp; ds.Bp', ds.C ];
          G = zeros(size(ds.B,2), 1);
          Ft = ds.Q' * ( (ds.Fp + ds.Fn) - ds.Ap * ds.U0');
          ds.F = [Ft;G];

        end

        function Kphi = getKphi(ds)
          for i = 1:ds.dof
              phi = zeros(ds.dof,1);
              phi(i) = 1;
              Aphi = gf_asm( 'volumic',...
                    [ 'm=data$1(#2);',...
                      't=comp(vGrad(#1).vGrad(#1).Base(#2));',...
                      'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
                      ' +t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
                      'M$1(#1,#1)+=sym(e(:,i,j,:,i,j,k).m(k));' ], ...
                      ds.mim, ds.mfu, ds.mfd, phi );
              %Kphi{i} = [Aphi, ds.Bp; ds.Bp', ds.C];
              Kphi{i} = Aphi;
          end
        end

        function L = getAdjointStiffness( ds, U )

          L = gf_asm('volumic',[ ...
                     'u=data(#2);',...
                     't=comp(vGrad(#2).vGrad(#2).Base(#1));',...
                     'e=(t{:,2,3,:,5,6,:}+t{:,3,2,:,5,6,:}',...
                     '+t{:,2,3,:,6,5,:}+t{:,3,2,:,6,5,:})/2;',...
                     'g=e(:,i,j,:,i,j,:);M(#2,#1)+=g(j,:,:).u(j)'],...
                     ds.mim, ds.mfd, ds.mfu, U);

        end

    end
end


