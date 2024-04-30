clear all;

load HANK_SS grid_a_NS pf_c pf_n pf_sav dist;

figure;
subplot(311);
plot(grid_a_NS,dist(:,1),'b-','LineWidth',2);
xlim([grid_a_NS(1) 1.0]); %grid_a_NS(end)]);
title("生産性: l_{L}","FontWeight","Normal");
%ylabel("密度","FontWeight","Normal")
box on;
grid on;
set(gca,'Fontsize',12)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(312);
plot(grid_a_NS,dist(:,2),'b-','LineWidth',2);
xlim([grid_a_NS(1) 5.0]); %grid_a_NS(end)]);
title("生産性: l_{M}","FontWeight","Normal");
%ylabel("密度","FontWeight","Normal")
box on;
grid on;
set(gca,'Fontsize',12)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(313);
plot(grid_a_NS,dist(:,3),'b-','LineWidth',2);
xlim([grid_a_NS(1) 50]); %grid_a_NS(end)]);
title("生産性: l_{H}","FontWeight","Normal");
xlabel("資産","FontWeight","Normal")
%ylabel("密度","FontWeight","Normal")
box on;
grid on;
set(gca,'Fontsize',12)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')

print -depsc2 dist.eps

