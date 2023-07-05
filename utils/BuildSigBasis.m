function A = BuildSigBasis(D,r,b,t)
Nd = length(D);
Nr = length(r);
Nb = length(b);
Nt = length(t);
A = zeros([Nb,Nt,Nd,Nr]);
for i = 1:Nd
    for j = 1:Nr
        Di = D(i);
        Rj = r(j);
        S_ij = exp(-Di*b(:))*exp(-Rj*t(:)');
        A(:,:,i,j) = S_ij;
    end
end
A = reshape(A,Nb*Nt,Nd*Nr);
end