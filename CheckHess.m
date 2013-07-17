function CheckHess(fn,y,N)

% CheckHess(fn,y,N)
%
%   This function checks the analytic Hessian against
%   a forward difference, printing out the step length,
%   error, and the exponent p in the model
%          error=O(h^p).
%   The error should be O(h^2), i.e. p=2.
%
%   The optional input N is the number of steps; the
%   default value is 4.
%
%   fn can be a function or an fstruct.  See "help fstruct"
%   for details.
%
%   The forward difference is computed at the point y,
%   in a random direction.

if nargin<3
   N=4;
end

% Get a random direction:

dy=randn(size(y));
dy=(max(1,norm(y))/(100*norm(dy)))*dy;

% Get the function value, gradient, and Hessian at y:

[fy,gy,Hy]=feval(fn,y);

% Compute the directional derivative:

if isstruct(Hy)
   der=fseval(Hy,dy);
else
   der=Hy*dy;
end

% Test the forward difference:

e0=0;
fprintf('h         error      p:\n')
disp('---------------------------------')
for h=2.^(0:-1:1-N)

   [f1,g1]=feval(fn,y+h*dy);
   e=norm(g1-gy-h*der);
   if h<1
      p=log(e0/e)/log(2);
      fprintf( '%.4e %.4e %.4e\n',h,e,p)
   else
      fprintf( '%.4e %.4e\n',h,e)
   end
   e0=e;

end
