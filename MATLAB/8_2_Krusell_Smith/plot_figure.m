% Purpose:
%  Plot result.
% 
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 09/11/2019  T. Yamada   Original code

clear;
clear global;
close all;
format short;

load results.mat;

global Params

%% CHANGE FONT FOR AXIS, TEXT AND LEGEND
set(0, 'defaultAxesFontSize', 16);
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontSize', 16);
set(0, 'defaultTextFontName', 'Arial');

%% MAKE FIGURE: POLICY FUNCTION

% policy function over current asset
figure;
plot(Params.grid,policy(:,1), '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(Params.grid,policy(:,2), '--', 'color', 'red', 'linewidth', 3); hold('off');
grid on;
title('Policy Function', 'fontsize', 16, 'fontweight', 'light');
xlabel('Current Asset', 'Fontsize', 16);
ylabel('Next Asset', 'Fontsize', 16);
legend('employed', 'unemployed', 'Location', 'NorthWest');
saveas (gcf, 'policy.eps', 'epsc2');


%% MAKE FIGURE: INVARIANT DISTRIBUTION

figure;
plot(Params.grid,dist, 'linewidth', 5);
grid on;
title('Distribution Function', 'fontsize', 16, 'fontweight', 'light')
xlabel('Asset')
ylabel('Density')
saveas (gcf,'distribution,eps','epsc2')

return;
