function [X,f] = MERDD_ProjectS(S,H_mu2sig,eps)
% check if the hankel forms of each row and column of M are positive
% definite or not, if not find the closet one 
[nb,nt] = size(S);
if(~exist('eps','var'))
    eps = 0;
end
M = H_mu2sig\S(:);
M = reshape(M,size(S));
min_eig = CheckMinEig(M);
if(min_eig>eps)
    X = S;
    f = 1;
else
    f = 0;
    cvx_begin quiet
        variable X(nb,nt)
        minimize(sum((H_mu2sig*X(:)-S(:)).^2))
        subject to 
        for i = 1:nt
            if(mod(nb,2))
                hankel(X(1:(nb+1)/2,i),X((nb+1)/2:end,i)') -eps*eye((nb+1)/2) == semidefinite((nb+1)/2);
                hankel(X(2:(nb+1)/2,i)-X(3:(nb+3)/2,i),X((nb+1)/2:end-1,i)'-X((nb+3)/2:end,i)') -eps*eye((nb-1)/2) == semidefinite((nb-1)/2);
            else
                hankel(X(2:(nb/2+1),i),X((nb/2+1):end,i)') -eps*eye(nb/2)== semidefinite(nb/2);
                hankel(X(1:nb/2,i)-X(2:(nb/2+1),i),X(nb/2:end-1,i)'-X((nb/2+1):end,i)') -eps*eye(nb/2)== semidefinite(nb/2);
            end
        end
        for j = 1:nb
            if(mod(nt,2))
                hankel(X(j,1:(nt+1)/2)',X(j,(nt+1)/2:end)) -eps*eye((nt+1)/2)== semidefinite((nt+1)/2);
                hankel(X(j,2:(nt+1)/2)'-X(j,3:(nt+3)/2)',X(j,(nt+1)/2:end-1)-X(j,(nt+3)/2:end)) -eps*eye((nt-1)/2)== semidefinite((nt-1)/2);
            else
                hankel(X(j,2:(nt/2+1))',X(j,(nt/2+1):end)) -eps*eye(nt/2-1) == semidefinite(nt/2-1);
                hankel(X(j,1:nt/2)'-X(j,2:(nt/2+1))',X(j,nt/2:end-1)-X(j,(nt/2+1):end)) -eps*eye(nt/2)== semidefinite(nt/2);
            end
        end

    cvx_end
    X = H_mu2sig*X(:);
    X = reshape(X,size(S)); 
end


end



function min_eig = CheckMinEig(X)
min_eig = 100;
[nb,nt]= size(X);
for i = 1:nt
    if(mod(nb,2))
        min_eig = min([min_eig, min(eig(hankel(X(1:(nb+1)/2,i),X((nb+1)/2:end,i)')))]);
        min_eig = min([min_eig, min(eig(hankel(X(2:(nb+1)/2,i)-X(3:(nb+3)/2,i),X((nb+1)/2:end-1,i)'-X((nb+3)/2:end,i)')))]); 
    else
        min_eig = min([min_eig, min(eig(hankel(X(2:(nb/2+1),i),X((nb/2+1):end,i)')))]);
        min_eig = min([min_eig, min(eig(hankel(X(1:nb/2,i)-X(2:(nb/2+1),i),X(nb/2:end-1,i)'-X((nb/2+1):end,i)')))]);
    end
end
for j = 1:nb
    if(mod(nt,2))
        min_eig = min([min_eig, min(eig(hankel(X(j,1:(nt+1)/2)',X(j,(nt+1)/2:end))))]);
        min_eig = min([min_eig, min(eig(hankel(X(j,2:(nt+1)/2)'-X(j,3:(nt+3)/2)',X(j,(nt+1)/2:end-1)-X(j,(nt+3)/2:end))))]);
    else
        min_eig = min([min_eig, min(eig(hankel(X(j,2:(nt/2+1))',X(j,(nt/2+1):end))))]);
        min_eig = min([min_eig, min(eig(hankel(X(j,1:nt/2)'-X(j,2:(nt/2+1))',X(j,nt/2:end-1)-X(j,(nt/2+1):end))))]);
    end
end
end