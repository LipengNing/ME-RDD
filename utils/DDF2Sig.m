function S = DDF2Sig(p,spec)
b = spec.b;
nb = length(b);
t = spec.t;
nt = length(t);
S = zeros(nb,nt);
D = spec.D;
R = spec.R;
delta_D = D(2)-D(1);
delta_R = R(2)-R(1);
for i = 1:nb
    for j = 1:nt
        X = exp(-b(i)*D(:))*exp(-t(j)*R(:)');
        S(i,j) = sum(X(:).*p(:))*delta_D*delta_R;
    end
end

end