function [p,D,R1,R2] = SpectralEst_r12dMUSIC(S,spec)
%multidimensional MUSIC algorithm based on block-hankel matrices
% for joint R1-T2 relaxation distribution 

[nb,nt2,nt1] = size(S); %
R1 = spec.R1;
R2 = spec.R2; %grid points for evaluating the spectra
D = spec.D;

N1 = length(R1);
N2 = length(R2);
ND = length(D);

ti = spec.ti; 
te = spec.te; te = te(:) - min(te);
b = spec.b; 

s1 = spec.s1; % number of columns for Hankel in t-space
s2 = spec.s2; % number of columns for Hankel in b-space
sb = spec.sb; %

ns = spec.ns; % expected number of signals

Hb = [];

if(~isfield(spec,'project'))
    spec.project = 1;%project=0: no projection; projection=1:only in T1 space; projection = 2 all space
end

if(~isfield(spec,'isdata'))
    spec.isdata=1;
end

if(spec.isdata)% the input may be a block Hankel matrix
    for i = 1:nb
        if(spec.project == 0)
            pp = 0;
        else
            pp = 1;
        end
        Sb = squeeze(S(i,:,:))';%t2xt1
        H = SpectralEst_Sig2r12Hankel(Sb,s1,s2,pp);    
        Hb = cat(3,Hb,H);%R12 Hankle for each b-value
    end
else
    Hb = S;
end

%%
hb = hankel(1:sb(1),sb(1):(sb(1)+sb(2)-1));

h_size = [s2(1)*s1(1)*sb(1)  s2(2)*s1(2)*sb(2)];%R1-space Hankel has an additional column

H = zeros(h_size);

for i = 1:nb
    Hi = double(hb == i);%indicator for non-zero hankel blocks
    H = H + kron(Hi,Hb(:,:,i));
end
%%%



nh = size(H,1);

[U,~,~]= svd(H);
U_n = U(:,(ns+1):end);


if(spec.project == 2)
    e = ones(nh,1);
    P = eye(nh) - e/(e'*e)*e';
else
    P = eye(nh);
end



p = zeros(ND,N2,N1);

ti_x = [0:s1(1)-1]'*(ti(2)-ti(1));
te_x = [0:s2(1)-1]'*(te(2)-te(1));
b_x = [0:sb(1)-1]'*(b(2)-b(1));

nb = length(b_x);
ne = length(te_x);
ni = length(ti_x);

if(spec.project == 0)% no projections
    PP = eye(ne*ni*nb);
else
    e = ones(ni,1);
    Pi = eye(ni) - e/(e'*e)*e';
    PP = kron(eye(nb),kron(eye(ne),Pi));
end


TE_x = kron(ones(nb,1),kron(te_x(:),ones(ni,1)));
B_x =  kron(b_x(:),ones(ne*ni,1));
TI_x = kron(ones(nb,1),kron(ones(ne,1),ti_x(:)));




DD = repmat(D,[1 N2*N1]);
RR2 = repmat(kron(R2,ones(1,ND)),[1 N1]);
RR1 = kron(R1,ones(1,ND*N2));
%tic
Basis = exp(-TE_x(:)*RR2-B_x(:)*DD).*(PP*exp(-TI_x(:)*RR1));

M = P*(U_n*U_n' + 1e-10*eye(nh))*P;

p = sum(Basis.*(P*Basis),1)./(sum(Basis.*(M*Basis),1));
p = reshape(p,ND,N2,N1);
%toc
% 
% tic
% parfor i = 1:N1
%     for j = 1:N2
%         for k = 1:ND
%             x = kron(exp(-b_x*D(k)),kron(exp(-te_x*R2(j)),Pi*exp(-ti_x*R1(i))));
%             p(k,i,j) = (x'*P*x)/(x'*P*(U_n*U_n' + 1e-10*eye(nh))*P*x);
%         end
%     end
% end
% toc


end



function H = SpectralEst_Sig2r12Hankel(S,s1,s2,project)
%transform each row of S into a Hankel matrix with k columns
%s1, s2: the number of rows (s1(1),s1(1)) and columns (s2(2), s2(2)) of Hankel in r1 and r2 spaces
%
%original was n_b x n_te
%S: (n_ti x n_te)
%sb(1)+sb(2) = nt+1;
if(nargin==3)
    project = 1;
end
[n1,n2] = size(S); %multi channel (e.g. TE data with different b-values)
%st(1)+st(2)==nt+1 
%sb(1)+sb(2)==nb+1
%

h_size = [s2(1)*s1(1)  s2(2)*s1(2)];%R1-space Hankel has an additional column

H1 = zeros(s1(1),s1(2),n2);%build a hankle in R1 for each column

e = ones(s1(1),1);


if(project)
    Pe = eye(s1(1)) - e/(e'*e)*e';
else
    Pe = eye(s1(1));
end

for i = 1:n2
    H1(:,:,i) =  Pe*hankel(S(1:s1(1),i),S((n1-s1(2)+1):n1,i));
end

hb = hankel(1:s2(1),s2(1):(s2(1)+s2(2)-1));

H = zeros(h_size);

for i = 1:n2
    Hi = double(hb == i);%indicator for non-zero hankel blocks
    H = H + kron(Hi,H1(:,:,i));
end


end