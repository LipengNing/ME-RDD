function [Mu, M_mu2sig]=MERDD_Sig2Mu_2d(S,delta_b,delta_t,D0,r0)
% Mu=MaxEntDDF_Sig2Mu_2d(S,delta_b)
% convert DWI signals to the moments of f(x) function
% input: S: dwi signal measured at b=0, delta_b, 2*delta_b, ...
% D0: upper bound of the diffusivity (default D_0= 3 mu^2/ms)
% output: Mu: the normalized moment
% M_mu2sig: transform matrix from mu to sig 

alpha = exp(-delta_b*D0);
beta = exp(-delta_t*r0);
[Nb,Nt] = size(S);
M_mu2sig = zeros(Nb,Nt,Nb,Nt);
sc = (1-alpha)*(1-beta)/(delta_b*delta_t);

for k = 1:Nb %k 0 to N-1
    for l = 1:Nt
        M = zeros(Nb,Nt);
        for i = 1:k
            for j = 1:l
                c1 = nchoosek(k-1,i-1)*alpha^(k-i)*(1-alpha)^(i-1);
                c2 = nchoosek(l-1,j-1)*beta^(l-j)*(1-beta)^(j-1);
                M(i,j) = c1*c2;
            end
        end
        M_mu2sig(k,l,:,:) = M;
    end
end
M_mu2sig = reshape(M_mu2sig,Nb*Nt,Nb*Nt);
M_mu2sig = sc*M_mu2sig;
Mu = M_mu2sig\S(:);
Mu = reshape(Mu,Nb,Nt);
        
