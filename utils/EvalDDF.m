function [peak_cor,center_cor,vol_frac] = EvalDDF(p,D,R,L)
%find the peak corrdinate and the vol_fraction for a DDF function with a
%given label map L (e.g. obtained using watershed of the ground truth)

id = unique(L(L>0));
K = length(id);
peak_cor = zeros(K,2);
center_cor = zeros(K,2);% center of mass
f = zeros(K,1);
for i = 1:K
    pi = zeros(size(p));
    pi(L==id(i)) = p(L==id(i));
    [~,ind] = max(pi(:));
    [x,y] = ind2sub(size(pi),ind);
    peak_cor(i,:) = [D(x) R(y)];
    f(i) = sum(pi(:));
%   center of mass
    ID_i = find(L==id(i));
    pii = p(ID_i);
    [X,Y] = ind2sub(size(pi),ID_i);
    Dx = D(X);
    Ry = R(Y);
    center_d = sum(pii(:).*Dx(:))/sum(pii);
    center_r = sum(pii(:).*Ry(:))/sum(pii);
    center_cor(i,:) = [center_d center_r];
end

vol_frac = f/sum(f(:))*100;
end