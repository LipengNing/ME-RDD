function ShowDDF(p,D,r)
figure;
delta_D = D(2)-D(1);
delta_r = r(2) - r(1);
pp = p/sum(p(:))/(delta_D*delta_r);
pp = pp(D<=2,r<=60);
D = D(D<=2);
r = r(r<=60);
Nd = length(D);
Nr = length(r);
%imagesc(r,D,pp);
h=surf(r,D,pp,'LineStyle',':');
%contour(ones(Nd,1)*r,D(:)*ones(1,Nr),pp);
set(gca,'CameraPosition',[250 7 3],'CameraViewAngle',10);
pD = sum(pp,2)*0.01;
pr = sum(pp,1)*0.01;

hold on;plot3(r,zeros(1,Nr),pr,'linewidth',2.5);
plot3(zeros(1,Nd),D,pD,'linewidth',2.5);
xlim([0 60]);
ylim([0 2]);

zlim([0 1.01*max([max(pp(:)) max(pr(:)) max(pD(:))])]);
set(gca,'fontsize',20);
set(gcf,'color','white');
zticklabels([]);
xlabel('$R~ (ms^{-1})$','interpreter','latex','Rotation',35);

ylabel('$D~ (\mu m^2/ms)$','interpreter','latex','Rotation',-35);


end