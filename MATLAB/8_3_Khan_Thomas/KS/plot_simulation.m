clear all;

load lumpyeqksresult.mat
load lumpyeqksresult2.mat

Z = [0.9328 0.9658 1.0000 1.0354 1.0720]';
GAMY  = 1.0160;
BETA  = 0.9540;
DELTA = 0.06;

figure;
subplot(313)
plot(Z(izvec),'r-','LineWidth',2.0)
xlim([1 2500]);
ylabel('A', 'FontSize',14);
xlabel('期間', 'FontSize',14);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(312)
plot(Kpvec,'b-','LineWidth',2.0)
xlim([1 2500]);
ylabel('K''', 'FontSize',14);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(311)
plot(Cvec,'b-','LineWidth',2.0)
xlim([1 2500]);
ylabel('C', 'FontSize',14);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
print -depsc2 KTsim.eps

Yvec = Cvec + Kpvec - (1-DELTA)*Kvec;
Xvec = Kpvec - (1-DELTA)*Kvec;
Rvec = [GAMY/BETA; Cvec(2:end)./Cvec(1:end-1)*GAMY/BETA]; % from BETA*C/C(+1)*R = 1

y = hpfilter(log(Yvec(501:end)),100);
c = hpfilter(log(Cvec(501:end)),100);
x = hpfilter(log(Xvec(501:end)),100);
n = hpfilter(log(Nvec(501:end)),100);
w = hpfilter(log(Wvec(501:end)),100);
r = hpfilter(log(Rvec(501:end)),100);

% disp([std(y)*100 std(c)/std(y) std(x)/std(y) std(n)/std(y) std(w)/std(y)]);
% disp([corr(y,c) corr(y,x) corr(y,n) corr(y,w)]);
disp([std(y)*100 std(c)/std(y) std(x)/std(y) std(n)/std(y) std(w)/std(y) std(r)/std(y)]);
disp([corr(y,c) corr(y,x) corr(y,n) corr(y,w) corr(y,r)]);