function ShowDDF2_invivo(p,D,r,D_max,r_max,trim,pmax)
if(~exist('D_max','var'))
    D_max=2;
end

if(~exist('r_max','var'))
    r_max = 60;
end

if(~exist('trim','var'))
    trim = 0;
end

if(isempty(trim))
    trim = 0;
end

if(~exist('pmax','var'))
    pmax = 0.04;
end


figure;
delta_D = D(2)-D(1);
delta_r = r(2) - r(1);

%pp = p/sum(p(:))/(delta_D*delta_r);
p = p(D<=D_max,r<=r_max);
D = D(D<=D_max);
r = r(r<=r_max);

pp = p/max(p(:))*1;



Nd = length(D);
Nr = length(r);
ppp = permute(pp,[2 1]);

%imagesc(D((1+trim):(end-trim)),r((1+trim):(end-trim)),ppp((1+trim):(end-trim),(1+trim):(end-trim)),[0 pmax]);
imagesc(D,r,ppp,[0 pmax]);

%imagesc(D(3:end-2),r(3:end-2),ppp(3:end-2,3:end-2),[0 .04]);
%h=surf(r(2:end),D(2:end),pp(2:end,2:end),'LineStyle',':');

%set(h,'LineWidth',.01);
%h = mesh(r,D,.3*pp);colormap(parula(10))

%h=contourf(ones(Nd,1)*r,D(:)*ones(1,Nr),pp,10);
%set(h,'linewidth',.1);
%mesh(r,D,pp);
set(gca,'CameraPosition',[14.417 593.369 1.312]);%,'CameraViewAngle',10);

pD = sum(pp,2);
pD = pD(3:(end-trim));
pr = sum(pp,1);
pr = pr(3:(end-trim));

pD = pD/sum(pD(:))*3;
pr = pr/sum(pr(:))*3;

m = max([max(pD(:)) max(pr(:))]);
if(m>0.2)
    pD = pD/m*0.19;
    pr = pr/m*0.19;
end
% m = max([max(pD(:)) max(pr(:))]);
% pD = pD/m*0.2;
% pr = pr/m*0.2;


hold on;plot3(D(3:(end-trim)),zeros(1,Nd-2-trim),pD*1.8,'linewidth',3);%3 10
plot3(zeros(1,Nr-2-trim),r(3:(end-trim)),pr*1.8,'linewidth',3);
xlim([0 D_max]);
ylim([0 r_max]);

%zlim([0 max([max(pr(:)) max(pD(:))])]);

%zlim([0 .2]);
zlim([0 0.2]);
set(gca,'fontsize',18,'XGrid','on','YGrid','on','ZGrid','on','Box','off');
set(gcf,'color','white');
zticklabels([]);
ylabel('$R~ (s^{-1})$','interpreter','latex','Rotation',30);

xlabel('$D~ (\mu m^2/ms)$','interpreter','latex','Rotation',-25);
set(gcf, 'Position', [0 0 770 700]);

end