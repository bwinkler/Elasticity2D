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
        function is = InverseSolver(ds, A0, Z, eps)
            is.ds = ds;
            is.A0 = A0;
            is.Z  = Z;
            is.eps = eps;
            %is.obj = OLSObjective( is.ds, is.Z, is.eps );
            is.obj = MOLSObjective( is.ds, is.Z, is.eps );
        end

        function [As, hist, cost, Ahist] = solve(is)
            f = @(A) is.obj.evaluate(A);

            CheckGrad( f, is.A0, 10 );

            [As, hist, cost, Ahist] = cgtrust1(is.A0, f, [is.tol, 0.1, 500, 500]);

            %[As, hist, cost] = gradproj(is.A0, f, 0.4 * ones(is.ds.dof,1),...
            %                            0.1* ones(is.ds.dof, 1), is.tol, 1E5);
            %Ahist = 0;
            
            %p.method = 1; 

            %p.alpha=0.1;
            %p.beta=0.8;

            %p.gamma = 0.9;
            %p.xsi = 0.8;
            %p.amin = 1E-4;

            %p.eps = 1E-10;

            %[As,hist] = marcotte( is.A0,...
            %                        f,...
            %                        2.6*ones(is.ds.dof,1),...
            %                        0.5*ones(is.ds.dof,1),...
            %                        is.tol, is.tol, ...
            %                        1E4,...
            %                        p);
            %cost = 0;
        end
    end

end

    
