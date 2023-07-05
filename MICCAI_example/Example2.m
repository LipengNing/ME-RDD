clear all
close all

%% subcortical
addpath ../
addpath ../utils


load InvivoData.mat
% configurate parameters
spec2d.maxIter=3000;
spec2d.r0 = 120;
spec2d.D0 = 3;
[~,~,nb,nt] = size(S);

spec2d.delta_t = 0.03;
spec2d.delta_b = 0.7;
spec2d.Nd = 500;
spec2d.Nr = 500;
spec2d.Tmin = 0.071;

spec2d.D = linspace(0,spec2d.D0,spec2d.Nd);
spec2d.R = linspace(0,spec2d.r0,spec2d.Nr);
spec2d.b = [0:nb-1]*spec2d.delta_b;
spec2d.t = spec2d.Tmin+[0:nt-1]*spec2d.delta_t;

spec2d.D = linspace(0,spec2d.D0,spec2d.Nd);
spec2d.R = linspace(0,spec2d.r0,spec2d.Nr);
spec2d.b = [0:nb-1]*spec2d.delta_b;
spec2d.t = spec2d.Tmin+[0:nt-1]*spec2d.delta_t;

%% extract signals from three voxels of different brain regions

S_gm  = squeeze(S(40,41,:,:));
S_wm  = squeeze(S(34,28,:,:));
S_subcort = squeeze(S(42,20,:,:));

%% GM

[p1,lambda1,D,r]=MaxEntDDF_Sig2DDF_2d_Delta1(S_gm,spec2d);
ShowDDF2_invivo(p1,D,r,3,120,2);

[p2,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta2(S_gm,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);

[p3,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(S_gm,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);

%% WM

[p1,lambda1,D,r]=MaxEntDDF_Sig2DDF_2d_Delta1(S_wm,spec2d);
ShowDDF2_invivo(p1,D,r,3,120,2);

[p2,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta2(S_wm,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);

[p3,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(S_wm,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);

%% subcortical GM

[p1,lambda1,D,r]=MaxEntDDF_Sig2DDF_2d_Delta1(S_subcort,spec2d);
ShowDDF2_invivo(p1,D,r,3,120,2);

[p2,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta2(S_subcort,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);

[p3,hatp2,lambda2,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(S_subcort,spec2d);
ShowDDF2_invivo(p2,D,r,3,120,0,0.001);



