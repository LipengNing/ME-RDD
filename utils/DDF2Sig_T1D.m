function S = DDF2Sig_T1D(p,spec)
%long TR and fixed TE
b = spec.b;
nb = length(b);
ti = spec.ti;
nt = length(ti);
S = zeros(nb,nt);
D = spec.D;
R1 = spec.R1;
delta_D = D(2)-D(1);
delta_R = R1(2)-R1(1);

for i = 1:nb
    for j = 1:nt
        X = exp(-b(i)*D(:))*abs(1-2*exp(-t(j)*R1(:)'));
        S(i,j) = sum(X(:).*p(:))*delta_D*delta_R;
    end
end

end