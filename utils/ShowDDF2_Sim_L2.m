function ShowDDF2_Sim_L2(p,D,r,D_max,r_max)
if(~exist('D_max','var'))
    D_max=2;
end

if(~exist('r_max','var'))
    r_max = 70;
end

figure;
delta_D = D(2)-D(1);
delta_r = r(2) - r(1);

%pp = p/sum(p(:))/(delta_D*delta_r);
p = p(D<=D_max,r<=r_max);
D = D(D<=D_max);
r = r(r<=r_max);

pp = p/max(p(:))*0.15;

ppp = permute(pp,[2 1]);

Nd = length(D);
Nr = length(r);
imagesc(D,r,ppp,[0 0.08]);
%h=surf(r(2:end),D(2:end),pp(2:end,2:end),'LineStyle',':');
%h=surf(r,D,pp,'LineStyle',':');

%set(h,'LineWidth',.01);
%h = mesh(r,D,.3*pp);colormap(parula(10))

%h=contourf(ones(Nd,1)*r,D(:)*ones(1,Nr),pp,10);
%set(h,'linewidth',.1);
%mesh(r,D,pp);
%set(gca,'CameraPosition',[291.7 9.85 1.55],'CameraViewAngle',10);
set(gca,'CameraPosition',[10.355 361.516 1.22]);%,'CameraViewAngle',10);


pD = sum(pp,2);
pr = sum(pp,1);
pD = pD/sum(pD(:))*5;
pr = pr/sum(pr(:))*5;

m = max([max(pD(:)) max(pr(:))]);
if(m>0.2)
    pD = pD/m*0.19;
    pr = pr/m*0.19;
end
% m = max([max(pD(:)) max(pr(:))]);
% pD = pD/m*0.2;
% pr = pr/m*0.2;


hold on;plot3(zeros(1,Nr),r,pr*2.5,'linewidth',3);
plot3(D,zeros(1,Nd),pD*2.5,'linewidth',3);
ylim([0 r_max]);
xlim([0 D_max]);

%zlim([0 max([max(pr(:)) max(pD(:))])]);

zlim([0 .2]);

set(gca,'fontsize',18,'XGrid','on','YGrid','on','ZGrid','on','Box','off');
set(gcf,'color','white');
zticklabels([]);
ylabel('$R~ (s^{-1})$','interpreter','latex','Rotation',30);

xlabel('$D~ (\mu m^2/ms)$','interpreter','latex','Rotation',-30);
set(gcf, 'Position', [0 0 700 550]);

end