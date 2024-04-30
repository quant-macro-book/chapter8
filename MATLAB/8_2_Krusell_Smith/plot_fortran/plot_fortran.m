clear;
close all;
format long;

%% CHANGE FONT FOR AXIS, TEXT AND LEGEND 

set(0, 'defaultAxesFontSize', 16);
set(0, 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');

%% LOAD DATA FROM EXCEL FILE
[sim1, header1] = xlsread('KS_Results', 'UnempRate');
[policy, header2] = xlsread('KS_Results', 'Policy');
[sim2, header3] = xlsread('KS_Results', 'LawOfMotion');
[appagg, header4] = xlsread('KS_Results', 'AppAgg');

%% UNEMPLOYMENT RATE

figure;
plot(sim1(5000:5100, 1), sim1(5000:5100, 3).*100, '-', 'color', 'blue', 'linewidth',3);
grid on;
xlim([5000, 5100]);
ylim([0, 12]);
xlabel('期間');
ylabel('失業率 (%)');
saveas (gcf,'Fig8_sim_unemp.eps','epsc2');

%% POLICY FUNCTION

figure;
plot(policy(:, 5), policy(:, 1), '-', 'color', 'blue', 'linewidth',3); hold('on');
plot(policy(:, 5), policy(:, 2), '--', 'color', 'red', 'linewidth',3);
plot(policy(:, 5), policy(:, 5), ':', 'color', 'black', 'linewidth',1); hold('off');
grid on;
xlim([0, 10]);
ylim([0, 10]);
xlabel('今期の資産：a');
ylabel("次期の資産：a'");
legend('就業', '失業', '45度線', 'Location', 'SouthEast');
saveas (gcf,'Fig8_policy.eps','epsc2');

%% SIMULATION

figure;
plot(sim2(5000:5100, 1), sim2(5000:5100, 3), '-', 'color', 'blue', 'linewidth',3);
grid on;
xlim([5000, 5100]);
ylim([34.5, 37]);
xlabel('期間');
ylabel('総資本：K');
saveas (gcf,'Fig8_sim_capital.eps','epsc2');

tfp = zeros(11000, 1);
for i = 1:11000
    if sim2(i, 2) == 1
        tfp(i) = 1.01;
    else
        tfp(i) = 0.99;
    end
end

figure;
plot(sim2(5000:5100, 1), tfp(5000:5100, 1), '-', 'color', 'blue', 'linewidth',3);
grid on;
xlim([5000, 5100]);
ylim([0.95, 1.05]);
xlabel('期間');
ylabel('TFP');
saveas (gcf,'Fig8_sim_tfp.eps','epsc2');

%% APPROXIMATE AGGREGATION

figure;
plot(appagg(:, 1), appagg(:, 2), '-', 'color', 'blue', 'linewidth',3); hold('on');
plot(appagg(:, 1), appagg(:, 3), '--', 'color', 'red', 'linewidth',3); hold('off');
grid on;
xlim([30, 40]);
ylim([30, 40]);
xlabel('現在の総資本：K');
ylabel("次期の総資本：K'");
legend('好況', '不況', 'Location', 'SouthEast');
saveas (gcf,'Fig8_app_agg.eps','epsc2');

return
