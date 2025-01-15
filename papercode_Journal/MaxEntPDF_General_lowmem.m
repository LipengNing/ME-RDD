function [p,lambda,s_est]=MaxEntPDF_General_lowmem(S,spec)
% A general formulation for MaxEntPDF suitable for any-dimensional data
% with provided spec.Basis functions
% find lambda that minimizes the following energy fuction
% Delta: =  int_0^D0 exp(-sum_i=1^(N) lambda_i
% (i-1)*delta_b*D)+\sum_i=1^(N) lambda_i*S(i+1)


S = S(:);
L = length(S(:));
lambdak = zeros(L,1); %initial value
maxIter = spec.maxIter;
if(~isfield(spec,'Basis'))
    error("The spec.Basis matrix is not provided!")
end

dTheta = spec.dTheta;%dTheta can be a vector for nonuniform samples, e.g. logspace

% if(length(dTheta)==1)
%     dTheta = dTheta*ones(size(spec.Basis,2),1);
%     spec.dTheta = dTheta;
% end

deltak = eval_objfunc(lambdak,S,spec);
stop = 0;
iter = 0;
beta = 0.5; % armijo parameter
sig = 1e-3; % armijo parameter
while(~stop)
    
    [g,H] = eval_gradients(lambdak,S,spec);
    
    if(min(eig(H))<1e-4)
        H = H + 1e-4*eye(size(H));
    end
    d = -H\g; % newton direction

   armijo_shrink = 1;
    alpha = 1;
    while(armijo_shrink)
        lambdak1 = lambdak+alpha*d;
        deltak1 = eval_objfunc(lambdak1,S,spec);
        if(deltak-deltak1> -sig*alpha*g'*d)
            armijo_shrink=0;
        else
            alpha = alpha*beta;
            if(norm(alpha*d)<1e-4) %step size too small
                armijo_shrink = 0;
                stop = 1;
            else
                if(norm(g)<1e-7*norm(S))%fitting error < 0.01%
                    armijo_shrink = 0;
                    stop = 1;
                end
            end

        end
    end
    %[norm(g) iter eval_objfunc(lambdak1,S,spec)]
    %iter
    %eval_objfunc(lambdak1,S,spec)
    lambdak = lambdak1;
    deltak = deltak1;
    iter = iter+1;
    if(iter>maxIter)
        stop = 1;
    end
    
    
end

lambda = lambdak;
p = eval_ddf(lambda,spec);
s_est = spec.Basis*(p(:).*dTheta(:));

end




function delta = eval_objfunc(lambda,S,spec)


dTheta = spec.dTheta; % unit interval/area
mu = spec.mu;

first_term = sum(exp(-lambda'*spec.Basis)*dTheta(:));
second_term = lambda'*S;
third_term = mu*(lambda'*lambda)/2;
delta = first_term + second_term + third_term;
end


function p = eval_ddf(lambda,spec)

p= exp(-lambda(:)'*spec.Basis);
end


function [g,H] = eval_gradients(lambda,S,spec)

p = eval_ddf(lambda,spec);

mu = spec.mu;
dTheta = spec.dTheta;
N = length(lambda);

Basis_p = spec.Basis.*repmat((p(:)*dTheta)',[N,1]);

g = -sum(Basis_p,2) +S+mu*lambda; % gradient gk = E(exp(-k*delta_b*D))

H = spec.Basis*Basis_p';



% parfor i = 1:L
%     H = H + spec.Basis(:,i)*spec.Basis(:,i)'*p(i)*dTheta(i);
% end

H = H+mu*(eye(N));
end



        
