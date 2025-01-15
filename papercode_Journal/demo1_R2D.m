%two dimensional R2-D distribution

te = ones(6,1)*[0.07 0.1 0.13 0.16 0.19];
b = [0:1:5]'*ones(1,5);

SNR = 200;
% model parameters
x = [18  .7 1;
     13  .3 1];%R2, D, vol-frac

S = zeros(size(b));
for i = 1:size(x,1)
    S = S + x(i,3)*exp(-x(i,1)*te-x(i,2)*b);
end

noise_level = mean(abs(S(:)))/SNR;

noise = noise_level*randn(size(S));
Sn = S + noise;
Sn = Sn/norm(Sn(:))*1e4;


%%  MaxEnt 
ND = 1000; NR = 1000;
D_ME = linspace(0.001,3,ND);
dD = gradient(D_ME); dD = dD/sum(dD)*3; 
R_ME = linspace(3,100,NR);
dR = gradient(R_ME); dR = dR/sum(dR)*97;
DD_ME = repmat(D_ME,[1 NR]);
RR_ME = kron(R_ME,ones(1,ND));

spec_ME.D = D_ME;
spec_ME.R = R_ME;
spec_ME.Basis = exp(-te(:)*RR_ME-b(:)*DD_ME);
spec_ME.dTheta = kron(dR,dD);
spec_ME.maxIter = 2000;
spec_ME.mu = 1e-4;

[p,lambda,s_est]= MaxEntPDF_General(Sn,spec_ME);
p = reshape(p,ND,NR);
figure;imagesc(R_ME(1:400),D_ME(1:400),p(1:400,1:400));
xlabel('R_2 (1/s)');
ylabel('D (\mum^2/ms)');
title('MaxEnt');

%% LR-LS+
ND_LS = 200; NR_LS = 200;
D_LS = linspace(0.001,3,ND_LS);
R_LS = linspace(3,100,NR_LS);
DD_LS = repmat(D_LS,[1 NR_LS]);
RR_LS = kron(R_LS,ones(1,ND_LS));
conf_LS.Basis = exp(-te(:)*RR_LS-b(:)*DD_LS);
conf_LS.mu = 1e-4;
conf_LS.maxIter = 2000;

[x,s_est] = LS_Nonnegative(Sn(:),conf_LS);
p = reshape(x/sum(x(:)),ND_LS,NR_LS);
figure;imagesc(R_LS(1:80),D_LS(1:80),p(1:80,1:80));
xlabel('R_2 (1/s)');
ylabel('D (\mum^2/ms)');
title('LR-LS_+');

%% HR-LS+
ND_LS = 1000; NR_LS = 1000;
D_LS = linspace(0.001,3,ND_LS);
R_LS = linspace(3,100,NR_LS);
DD_LS = repmat(D_LS,[1 NR_LS]);
RR_LS = kron(R_LS,ones(1,ND_LS));
conf_LS.Basis = exp(-te(:)*RR_LS-b(:)*DD_LS);
conf_LS.mu = 1e-4;
conf_LS.maxIter = 2000;

[x,s_est] = LS_Nonnegative(Sn(:),conf_LS);
p = reshape(x/sum(x(:)),ND_LS,NR_LS);
figure;imagesc(R_LS(1:400),D_LS(1:400),p(1:400,1:400));
xlabel('R_2 (1/s)');
ylabel('D (\mum^2/ms)');
title('LS_+');



%% MUSIC
param_MUSIC.D = D_ME;
param_MUSIC.R = R_ME;
param_MUSIC.ns = 2;
param_MUSIC.t = te(1,:)';
param_MUSIC.b = b(:,1);
param_MUSIC.st = [3 3];
param_MUSIC.sb = [4 3];
param_MUSIC.isdata = 1;

[p,D,R] = SpectralEst_rdMUSIC(Sn,param_MUSIC);
p = reshape(p/sum(p(:)),ND,NR);
figure;imagesc(R_ME(1:400),D_ME(1:400),p(1:400,1:400));
xlabel('R_2 (1/s)');
ylabel('D (\mum^2/ms)');
title('MUSIC');
