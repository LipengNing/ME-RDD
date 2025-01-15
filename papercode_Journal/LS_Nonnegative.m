function [x,y_est] = LS_Nonnegative(y,spec)
%ADMM method for solving min_X .5||Ax-y||^2 + 0.5*c ||x||^2, x>=0

A = spec.Basis;
% x0 = lsqnonneg(A,y);
c1 = spec.mu;
% if(c1 == 0) 
%     x = x0;
%     y_est = A*x;
% else
maxIter = spec.maxIter;
c2 = 1;% auxilliary parameter

c = c1 + c2;
[m,n] = size(A);

M = inv(A*A'+c*eye(m)); % matrix used to solve the inverse problem

%initialize augment and multiplier
xh_k = zeros(n,1);%initial value based on lsqnonneg
u_k = zeros(n,1);
stop = 0;

eps_pri = 1e-8;
eps_dual = 1e-8;
iter = 1;
    while(~stop)
        %update x
        
        x_tmp = A'*y+c2*(xh_k-u_k);
        x_k1 = (x_tmp - A'*(M*A*x_tmp))/c;%(eye(n)-A'*M*A)*x_tmp/c;
        
        %update xh
        xh_k1 = x_k1 + u_k;
        xh_k1(xh_k1<0) = 0;
        % update u_k
        u_k1 = u_k + x_k1 - xh_k1;
        
        if(norm(xh_k1-xh_k)<eps_dual)
            stop = 1;
        end
        if(norm(x_k1 -xh_k)<eps_pri)
            stop = 1;
        end
        if(iter>maxIter)
            stop = 1;
        else 
            iter = iter + 1;
        end
        u_k = u_k1;
        x_k = x_k1;
        xh_k = xh_k1;
        
    end
    
    x = xh_k;
    y_est = A*x(:);
    
end
%end
