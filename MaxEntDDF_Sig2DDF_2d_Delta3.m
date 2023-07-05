function [p,hatp,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(S,spec)
% [p,D]=MaxEntDDF_Sig2DDF_2d_Delta3(S,delta_b,D0)
% find lambda that minimizes the following energy fuction
% Delta: =  int_0^D0 exp(-sum_i=1^(N) lambda_i
% (i-1)*delta_b*D)+\sum_i=1^(N) lambda_i*S(i+1)
% p = hapt*exp(Tmin/delta_t*r);

delta_b = spec.delta_b;
delta_t = spec.delta_t;
D0=spec.D0;
r0 = spec.r0;
Nd = spec.Nd;
Nr = spec.Nr;
Tmin = spec.Tmin;
maxIter = spec.maxIter;
D = linspace(0,D0,Nd);
r = linspace(0,r0,Nr);
[Nb,Nt] = size(S);
lambdak = zeros(Nb,Nt); %initial value
deltak = eval_objfunc(lambdak,S,delta_b,delta_t,D,r);
stop = 0;
iter = 0;
beta = 0.5; % armijo parameter
sig = 1e-3; % armijo parameter
alpha = 1; %initial stepsize
while(~stop)
    [g,H] = eval_gradients(lambdak,S,delta_b,delta_t,D,r);
    d = -(H+1e-3*eye(size(H)))\g; % newton direction
    d = reshape(d,Nb,Nt);
    armijo_shrink = 1;
    while(armijo_shrink)
        lambdak1 = lambdak+alpha*d;
        deltak1 = eval_objfunc(lambdak1,S,delta_b,delta_t,D,r);
        if(deltak-deltak1> -sig*alpha*g'*d(:))
            armijo_shrink=0;
        else
            alpha = alpha*beta;
            if(norm(alpha*d)<1e-4) %step size too small
                armijo_shrink = 0;
                stop = 1;
            else
                if(norm(g)<1e-4*norm(S))
                    armijo_shrink = 0;
                    stop = 1;
                end
            end
        end
    end
%     norm(g)
%     iter
    %eval_objfunc(lambdak1,S,delta_b,delta_t,D,r)
    lambdak = lambdak1;
    deltak = deltak1;
    iter = iter+1;
    if(iter>maxIter)
        stop = 1;
    end
end

lambda = lambdak;
hatp = eval_ddf(lambda,S,delta_b,delta_t,D,r);
%hatp(hatp<1e-5) = 1e-5;
p = hatp.*(ones(Nd,1)*exp(Tmin*r));
end




function delta = eval_objfunc(lambda,S,delta_b,delta_t,D,r)
[Kb,Kt] = size(S);
deltaD = D(2)-D(1);
deltar = r(2)-r(1);
b = [0:Kb-1]'*delta_b;
t = [0:Kt-1]'*delta_t;

first_term = 0;
for k = 1:Kb
    for l = 1:Kt
        Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
        first_term = first_term + Mkl;
    end
end


%sec_term = ones(length(D),1)*exp(5*Tmin*r);
sec_term = exp(-delta_b*D(:))*exp(-delta_t*r);
first_term = exp(first_term).*sec_term;

first_term = sum(first_term(:))*deltaD*deltar;

second_term = lambda(:)'*S(:);

delta = first_term + second_term;
end


function p = eval_ddf(lambda,S,delta_b,delta_t,D,r)
[Kb,Kt] = size(S);
b = [0:Kb-1]'*delta_b;
t = [0:Kt-1]'*delta_t;
first_term = 0;
for k = 1:Kb
    for l = 1:Kt
        Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
        first_term = first_term + Mkl;
    end
end

%sec_term = ones(length(D),1)*exp(5*Tmin*r);
sec_term = exp(-delta_b*D(:))*exp(-delta_t*r);
p= exp(first_term).*sec_term;
end


function [g,H] = eval_gradients(lambda,S,delta_b,delta_t,D,r)
p = eval_ddf(lambda,S,delta_b,delta_t,D,r);
[Kb,Kt] = size(S);
deltaD = D(2)-D(1);
deltar = r(2)-r(1);
b = [0:Kb-1]'*delta_b;
t = [0:Kt-1]'*delta_t;
g = zeros(Kb,Kt);

Nd = length(D);
Nr = length(r);
HH = zeros(Nd,Nr,Kb,Kt);
for k = 1:Kb
    for l = 1:Kt
        M = exp(-b(k)*D')*exp(-t(l)*r);
        g(k,l) = -sum(M(:).*p(:))*deltaD*deltar+S(k,l);
        HH(:,:,k,l) = M;
    end
end
g = g(:);
HH = reshape(HH,Nd*Nr,Kb*Kt);
HH = HH.*repmat(sqrt(p(:)),1,Kb*Kt);
H = HH'*HH*deltaD*deltar;
end



        
