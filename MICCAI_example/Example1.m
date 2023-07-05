clear all
close all
addpath ../utils
addpath ..
spec.w = [2 5 3]; % weight of different components
spec.center = [1.5 10; 0.5 40; 1.5 40];
spec.b = 0.5*[0:10];
spec.t = .05+[0:6]*0.025;

spec.delta_b =  spec.b(2)-spec.b(1);
spec.delta_t = spec.t(2)-spec.t(1);
spec.D0 = 3;
spec.r0 = 80;
spec.Nd = 400;
spec.Nr = 400;
spec.Tmin = min(spec.t);
spec.maxIter = 1000;

Sig(:,:,1) = diag([0.01 5]);
Sig(:,:,2) = diag([0.01 5]);
Sig(:,:,3) = diag([0.01 5]);
spec.sig = Sig;

spec.D = linspace(0,3,500);
spec.R = linspace(0,100,400);
spec.Tmin = min(spec.t);


p = CreatDDF(spec);
ShowDDF2(p,spec.D,spec.R);
% %%

S = DDF2Sig(p,spec);
S = S*5000;

[p1,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta1(S,spec);
ShowDDF2_Sim(p1,D,r);

% export_fig('figs_sim/conf2_Delta1.pdf');
% 
[p2,hatp,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta2(S,spec);
ShowDDF2_Sim(p2,D,r);
% export_fig('figs_sim/conf2_Delta2.pdf');
% 
[p3,hatp3,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(S,spec);
ShowDDF2_Sim(p3,D,r);
% export_fig('figs_sim/conf2_Delta3.pdf');
% 
% 
% 
snr =300;
sig = S./snr;
% 
Si = S+sig.*randn(size(S));
%apply three  ME-RDD methods
[p1,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta1(Si,spec);
ShowDDF2_Sim(p1,D,r);
% 
% export_fig('figs_sim/conf2_Delta1_snr300.pdf');

[p2,hatp,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta2(Si,spec);
ShowDDF2_Sim(p2,D,r);
% export_fig('figs_sim/conf2_Delta2_snr300.pdf');
% 

[p3,hatp3,lambda,D,r]=MaxEntDDF_Sig2DDF_2d_Delta3(Si,spec);
ShowDDF2_Sim(p3,D,r);
% export_fig('figs_sim/conf2_Delta3_snr300.pdf');