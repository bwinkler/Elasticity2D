function [a1,err,k,y]=secant(f,a0,a1,delta,epsilon,maxIter)

for k=1:maxIter    
    a2=a1-f(a1)*(a1-a0)/(f(a1)-f(a0));    
    err=abs(a2-a1);
    relerr=2*err/(abs(a2)+delta);
    a0=a1;
    a1=a2;
    y=f(a1);
    if (err<delta)|(relerr<delta)|(abs(y)<epsilon)
      break
    end
end
