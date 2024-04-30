clear all;

yy = importdata("IRF.dat")
y1HANK = yy(:,1);
y1RANK = yy(:,2);
y2HANK = yy(:,3);
y2RANK = yy(:,4);
y3HANK = yy(:,5);
y3RANK = yy(:,6);

figure;
subplot(311);
plot(y1HANK,'k-',"LineWidth",2);
hold on;
plot(y1RANK,'r--');
xlim([1,40]);
legend("HANK","RANK");
title("消費");
subplot(312);
plot(y2HANK,'k-',"LineWidth",2);
hold on;
plot(y2RANK,'r--');
xlim([1,40]);
title("インフレ率");
subplot(313);
plot(y3HANK,'k-',"LineWidth",2);
hold on;
plot(y3RANK,'r--');
xlim([1,40]);
title("政策金利");