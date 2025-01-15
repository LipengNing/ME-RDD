%two dimensional R2-D distribution

b = [0 1 2 3 4 5]'*[1 1 1 1 1];
te = [1 1 1 1 1 1]'*[0.07 0.1 0.13 0.16 0.19];
ti = [300 500 700 900 1100 1300]/1000;
te = te(:,:,ones(1,6));
b = b(:,:,ones(1,6));
ti = kron(ti,ones(6,5));
ti = reshape(ti,[6 5 6]);
%model parameters
x = [18  .7  0.7 1;
     13  .3  1.1  1;
     30 0.2  1.5 1];%R2, D, vol-frac

S = zeros(size(b));
for i = 1:size(x,1)
    S = S + x(i,3)*exp(-x(i,1)*te-x(i,2)*b).*(1-2*exp(-x(i,3)*ti));
end

SNR = 200;
noise_level = mean(abs(S(:)))/SNR;
noise = noise_level*randn(size(S));
Sn = S + noise;
Sn = Sn/norm(Sn(:))*1e4;


%%  MaxEnt 
ND = 200; NR2 = 200; NR1=200;
D_ME = linspace(0.001,3,ND);
dD = gradient(D_ME); dD = dD/sum(dD)*3; 
R2_ME = linspace(3,100,NR2);
dR2 = gradient(R2_ME); dR2 = dR2/sum(dR2)*97;
R1_ME = linspace(0.3,3,NR1);%this is R1 
dR1 = gradient(R1_ME); dR1 = dR1/sum(dR1)*2.7; 
DD = repmat(D_ME,[1 NR2*NR1]);
RR2 = repmat(kron(R2_ME,ones(1,ND)),[1 NR1]);
RR1 = kron(R1_ME,ones(1,ND*NR2));

spec_ME.D = D_ME;
spec_ME.R2 = R2_ME;
spec_ME.R1 = R1_ME; 
spec_ME.Basis = exp(-te(:)*RR2-b(:)*DD).*(1-2*exp(-ti(:)*RR1));
spec_ME.dTheta = dR1(1)*dR2(1)*dD(1);
spec_ME.maxIter = 1000;
spec_ME.mu = 1e-4;
% MaxEntPDF_General_lowmem: low-memory version of MaxEntPDF_General
tic
[p,lambda,s_est]= MaxEntPDF_General(Sn,spec_ME);
t_me=toc
p = reshape(p,ND,NR2,NR1);
p = permute(p,[1 3 2]);
p_DR2 = squeeze(mean(p,2));
p_DR2=p_DR2/prctile(p_DR2(:),99.9);

p_DR1 = squeeze(mean(p,3));
p_DR1 = p_DR1/prctile(p_DR1(:),99.9);

p_R1R2 = squeeze(mean(p,1))';
p_R1R2 = p_R1R2/prctile(p_R1R2(:),99.9);

trueParam = x;

    
figure;
imagesc(R2_ME,D_ME,p_DR2,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MaxEnt');
set(gca,'fontsize',20);
xlim([5 45]);
ylim([0.001 1.2]);
ylabel('D (\mum^2/ms)');
xlabel('R_2 (1/s)');
set(gcf,'color','white');

figure;
imagesc(R1_ME,D_ME,p_DR1,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,3),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MaxEnt');
set(gca,'fontsize',20);
xlim([0.3 2]);
ylim([0.001 1.2]);
ylabel('D (\mum^2/ms)');
xlabel('R_1 (1/s)');
set(gcf,'color','white');


figure;
imagesc(R2_ME,R1_ME,p_R1R2',[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,3),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MaxEnt');
set(gca,'fontsize',20);
xlim([5 45]);
ylim([0.3 2]);
xlabel('R_2 (1/s)');
ylabel('R_1 (1/s)');
set(gcf,'color','white');


%% LS+

tic
[x,s_est] = LS_Nonnegative(Sn(:),spec_ME);
t_ls = toc
p = reshape(x/sum(x(:)),ND,NR2,NR1);

p = permute(p,[1 3 2]);

p_DR2 = squeeze(mean(p,2));
p_DR2=p_DR2/prctile(p_DR2(:),99.9);

p_DR1 = squeeze(mean(p,3));
p_DR1 = p_DR1/prctile(p_DR1(:),99.9);

p_R1R2 = squeeze(mean(p,1))';
p_R1R2 = p_R1R2/prctile(p_R1R2(:),99.9);

    
figure;
imagesc(R2_ME,D_ME,p_DR2,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('LS_+');
set(gca,'fontsize',20);
xlim([5 45]);
ylim([0.001 1.2]);
ylabel('D (\mum^2/ms)');
xlabel('R_2 (1/s)');
set(gcf,'color','white');

figure;
imagesc(R1_ME,D_ME,p_DR1,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,3),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('LS_+');
set(gca,'fontsize',20);
xlim([0.3 2]);
ylim([0.001 1.2]);
ylabel('D (\mum^2/ms)');
xlabel('R_1 (1/s)');
set(gcf,'color','white');


figure;
imagesc(R2_ME,R1_ME,p_R1R2',[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,3),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('LS_+');
set(gca,'fontsize',20);
xlim([5 45]);
ylim([0.3 2]);
xlabel('R_2 (1/s)');
ylabel('R_1 (1/s)');
set(gcf,'color','white');



%% MUSIC

param_MUSIC.D = D_ME;
param_MUSIC.R2 = R2_ME;
param_MUSIC.R1 = R1_ME;
param_MUSIC.ns = 3;
param_MUSIC.te = te(1,:)';
param_MUSIC.b = b(:,1);
param_MUSIC.ti = squeeze(ti(1,1,:));
param_MUSIC.s2 = [3 3];% size of hankle in t2
param_MUSIC.sb = [4 3];% size of hankle in diffusion
param_MUSIC.s1 = [4 3];% xx in t1
param_MUSIC.ns = 3;
param_MUSIC.project = 1;%project to null(ones) in the R1 space
tic
[p,D,R2,R1] = SpectralEst_r12dMUSIC(Sn,param_MUSIC);
t_music = toc
p = reshape(p/sum(p(:)),ND,NR2,NR1);
p = permute(p,[1 3 2]);

p_DR2 = squeeze(mean(p,2));
p_DR2=p_DR2/prctile(p_DR2(:),99.9);

p_DR1 = squeeze(mean(p,3));
p_DR1 = p_DR1/prctile(p_DR1(:),99.9);

p_R1R2 = squeeze(mean(p,1))';
p_R1R2 = p_R1R2/prctile(p_R1R2(:),99.9);



figure;
imagesc(R2_ME,D_ME,p_DR2,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MUSIC');
set(gca,'fontsize',20);
xlim([3 60]);
ylim([0.001 2]);
ylabel('D (\mum^2/ms)');
xlabel('R_2 (1/s)');
set(gcf,'color','white');


figure;
imagesc(R1_ME,D_ME,p_DR1,[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,3),trueParam(m,2),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MUSIC');
set(gca,'fontsize',20);
xlim([0.3 2]);
ylim([0.001 2]);
ylabel('D (\mum^2/ms)');
xlabel('R_1 (1/s)');
set(gcf,'color','white');


figure;
imagesc(R2_ME,R1_ME,p_R1R2',[0 1]);hold on;
for m = 1:size(trueParam,1)
    plot(trueParam(m,1),trueParam(m,3),'r+','MarkerFaceColor','r','MarkerSize',25,'LineWidth',1.5);
end
title('MUSIC');
set(gca,'fontsize',20);
xlim([3 60]);
ylim([0.3 2]);
xlabel('R_2 (1/s)');
ylabel('R_1 (1/s)');
set(gcf,'color','white');
