function [Sp,E] = MERDD_ProjSig(S,A)
%  [Sp,E] = MaxEntDDF_ProjSig(S,A)
% ADMM algorithm to project signal S to the convex set A*S>0 (eps)

[L,N] = size(A);
b=1e-3*ones(L,1);
MaxIter = 2000;
rho = 10;
M = (eye(N)+rho*A'*A)\eye(N);
stop = 0;
iter = 0;

x0 = S(:);
xk = x0;
yk = A*xk-b;
yk(yk<0) = 0;
lambdak = zeros(L,1);
while(~stop)
   %updata X
   xk1 = M*(x0+rho*A'*(b+yk)-A'*lambdak);
   %update yk
   yk1 = A*xk1-b+lambdak/rho;
   yk1(yk1<0) = 0;
   %update lambdak1
   lambdak1 = lambdak+rho*(A*xk1-b-yk1);
   %%check stopping
   res1 = norm(A*xk1-b-yk1);
   res2 = norm(xk1-xk);
   if(res1<1e-3&&res2<1e-3)
       xk = xk1;
       stop = 1;
   end
   iter = iter+1;
   xk = xk1;
   yk = yk1;
   lambdak = lambdak1;
   if(iter>MaxIter)
       stop = 1;
   end

end

Sp = reshape(xk,size(S));
E = S-Sp;

end