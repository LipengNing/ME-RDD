function Mu=MERDD_Sig2Mu(Data,spec)
% Mu=MaxEntDDF_Sig2Mu(Data,spec)
% convert DWI signals to the moments of f(x) function

mask = spec.mask;
delta_b = spec.delta_b;
D0 = spec.D0;
b = spec.bvals;
ub = unique(b);
nb = length(ub);
nu = sum(b==max(ub));% number of gradient directions

[nx,ny,nz,~] = size(Data);
Data = reshape(Data,nx*ny*nz,nb);
Data_b0 = mean(Data(:,b==0),2);
Data_b0 = Data_b0(:,ones(nu,1)); 
Data_nb0 = Data(:,b>0);
Data_nb0 = reshape(Data_nb0,nx*ny*nz,nu,nb-1);
Data = cat(3,Data_b0,Data_nb0);%nxnynz x nu x nb
Mu = zeros(size(Data));

for i = nx*ny*nz
    if(mask(i)>0)
        for j = 1:nu
            data = squeeze(Data(i,j,:));
            mu=MaxEntDDF_Sig2Mu_1d(data,delta_b,D0);
            Mu(i,j,:) = mu;
        end
    end
end
