function ShowDDF2(p,D,r,D_max,r_max)
if(~exist('D_max','var'))
    D_max=2;
end

if(~exist('r_max','var'))
    r_max = 60;
end

figure;
delta_D = D(2)-D(1);
delta_r = r(2) - r(1);

%pp = p/sum(p(:))/(delta_D*delta_r);
p = p(D<=D_max,r<=r_max);
D = D(D<=D_max);
r = r(r<=r_max);

pp = p/max(p(:))*0.15;



Nd = length(D);
Nr = length(r);
%imagesc(r,D,pp);
%h=surf(r(2:end),D(2:end),pp(2:end,2:end),'LineStyle',':');
h=surf(r,D,pp,'LineStyle',':');

set(h,'LineWidth',.01);
%h = mesh(r,D,.3*pp);colormap(parula(10))

%h=contourf(ones(Nd,1)*r,D(:)*ones(1,Nr),pp,10);
%set(h,'linewidth',.1);
%mesh(r,D,pp);
set(gca,'CameraPosition',[291.7 9.85 1.55],'CameraViewAngle',10);
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


hold on;plot3(r,zeros(1,Nr),pr,'linewidth',3);
plot3(zeros(1,Nd),D,pD,'linewidth',3);
xlim([0 r_max]);
ylim([0 D_max]);

%zlim([0 max([max(pr(:)) max(pD(:))])]);

zlim([0 .2]);

set(gca,'fontsize',18,'XGrid','on','YGrid','on','ZGrid','on','Box','off');
set(gcf,'color','white');
zticklabels([]);
xlabel('$R~ (ms^{-1})$','interpreter','latex','Rotation',35);

ylabel('$D~ (\mu m^2/ms)$','interpreter','latex','Rotation',-30);
set(gcf, 'Position', [0 0 1700 1000]);

end