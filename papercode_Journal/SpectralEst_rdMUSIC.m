function [p,D,R] = SpectralEst_rdMUSIC(S,spec)
%multidimensional MUSIC algorithm based on block-hankel matrices
% for, e.g., joint relaxation-diffusion MRI

[nb,nt] = size(S); %
D = spec.D;
R = spec.R; %grid points for evaluating the spectra
Nd = length(D);
Nr = length(R);

t = spec.t; t = t(:) - min(t);
b = spec.b; b = b(:) - min(b);

st = spec.st; % number of columns for Hankel in t-space
sb = spec.sb; % number of columns for Hankel in b-space

ns = spec.ns; % expected number of signals

if(spec.isdata)% the input may be a block Hankel matrix
    H = SpectralEst_Sig2rdHankel(S,st,sb);
else
    H = S;
end

nh = size(H,1);

[U,~,~]= svd(H);

U_n = U(:,(ns+1):end);
p = zeros(Nd,Nr);

b_x = [0:sb(1)-1]'*(b(2)-b(1));
t_x = [0:st(1)-1]'*(t(2)-t(1));


L = size(U_n,1);

if(~isfield(spec,'project'))
    spec.project = 0;
end

if(spec.project)
    e = ones(L,1);
    P = eye(L)- e/(e'*e)*e';
else
    P = eye(L);
end

nd = length(D);
nr = length(R);
DD = kron(ones(1,nr),D(:)');
RR = kron(R(:)',ones(1,nd));
nb = length(b_x);
nt = length(t_x);
B_x = kron(b_x(:),ones(nt,1));
TE_x = kron(ones(nb,1),t_x(:));

Basis = exp(-B_x*DD-TE_x*RR);
M = P*(U_n*U_n' + 1e-10*eye(nh))*P';
p = sum(Basis.*(P*Basis),1)./sum(Basis.*(M*Basis),1);
p = reshape(p,Nd,Nr);
% 
% for i = 1:Nd
%     for j = 1:Nr
%         x = kron(exp(-b_x*D(i)),exp(-t_x*R(j)));
%         p(i,j) = (x'*P*x)/(x'*P*(U_n*U_n' + 1e-10*eye(nh))*P*x);
%     end
% end

end


function H = SpectralEst_Sig2rdHankel(S,st,sb)
%transform each row of S into a Hankel matrix with k columns
%st, sb: the number of rows (st(1),sb(1)) and columns (st(2), sb(2)) of Hankel in t and b spaces
%

%sb(1)+sb(2) = nt+1;
[nb,nt] = size(S); %multi channel (e.g. TE data with different b-values)
%st(1)+st(2)==nt+1 
%sb(1)+sb(2)==nb+1
%

h_size = [sb(1)*st(1)  sb(2)*st(2)];

Ht = zeros(st(1),st(2),nb);
for i = 1:nb
    Ht(:,:,i) = hankel(S(i,1:st(1))',S(i,(nt-st(2)+1):nt));
end

hb = hankel(1:sb(1),sb(1):(sb(1)+sb(2)-1));

H = zeros(h_size);

for i = 1:nb
    Hi = double(hb == i);%indicator for non-zero hankel blocks
    H = H + kron(Hi,Ht(:,:,i));
end

end