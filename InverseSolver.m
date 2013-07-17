classdef InverseSolver < handle
    % Solves Stokes/Incompressibility inverse problem

    properties
        ds
        obj
        tol = 1E-10;
        A0
        Z
        eps = 1E-6;
    end
    methods
        function is = InverseSolver(ds, A0, Z, eps, objMeth, gradMeth)
            is.ds = ds;
            is.A0 = A0;
            is.Z  = Z;
            is.eps = eps;

            switch objMeth
              case 'MOLS'              
                is.obj = MOLSObjectiveJ2( is.ds, is.Z, is.eps, gradMeth );
              case 'MOLSJ2'              
                is.obj = MOLSObjectiveJ2( is.ds, is.Z, is.eps, gradMeth );
              case 'MOLSJ1'
                is.obj = MOLSObjectiveJ1( is.ds, is.Z, is.eps, gradMeth);
              case 'OLS'
                is.obj = OLSObjective( is.ds, is.Z, is.eps, gradMeth );
              case 'EE'
                 is.obj = EEObjective( is.ds, is.Z, is.eps);
              otherwise
                error('Invalid objective function method specified');
            end
        end

        function [As, hist, cost, Ahist] = solve(is)
            f = @(A) is.obj.evaluate(A);

            disp('Checking gradient...');
            CheckGrad(f, is.A0, 5);
            disp('Checking hessian...');
            CheckHess(f, is.A0, 10);

            [As, hist, cost, Ahist] = cgtrust1(is.A0, f, [is.tol, 0.1, 500, 500]);

            %[As, hist, cost, Ahist] = cgtrust1(is.A0, f, is.tol);

            %[As, hist, cost] = bfgswopt(is.A0, f, is.tol,100000);
            %Ahist = 0;
            % [As, hist, cost] = gradproj(is.A0, f, 3 * ones(is.ds.dof,1),...
            %                              1* ones(is.ds.dof, 1), is.tol, 1E5);
            %Ahist = 0;
        end
    end

end

    
