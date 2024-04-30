clear all;

% parameters
p.GAMMA = 2.0;
p.PSI = 2.0;
p.BL = 0.0;
p.MU = 1.2;
p.SUB_MON_COM = (p.MU-1.0)/p.MU;

p.TAU_BAR = [0 0 2.031292672482304]';
p.B4YSS = 5.5/4.0;

p.MAX_A = 75;
p.MIN_A = 0;
p.MAX_A_NS = 75;
p.MIN_A_NS = 0;
p.N_A = 200;
p.N_A_NS = 1000;
p.GRID_A_SHIFT = 0.354983310608304; % why this number?
p.BL = 0;

p.RHO_Z = 0.96566;
p.VAR_Z_INN = 0.01695^(1/2); % variance^(1/2)
p.N_Z = 3;

p.MAX_BETA_ITER = 50;
p.MAX_EGM_ITER = 1000;
p.MAX_EGM_CONST_ITER = 1000;
p.BETA_ERR_TOL = 1e-4; %1e-12;
p.EGM_ERR_TOL = 1e-8; %1e-15;
p.EGM_CONST_ERR_TOL = 1e-12;
p.NS_ERR_TOL = eps; %1e-18;

W_SS = 1.0/(p.MU*(1.0 - p.SUB_MON_COM)); % = 1
R_SS = 1.005;
PI_SS = 1.0;

% Discretize the endogenous state space for assets
grid_a_preshift = linspace(log(p.MIN_A_NS+p.GRID_A_SHIFT),log(p.MAX_A_NS+p.GRID_A_SHIFT),p.N_A)';
grid_a = exp(grid_a_preshift) - p.GRID_A_SHIFT;

grid_a_preshift_NS = linspace(log(p.MIN_A+p.GRID_A_SHIFT),log(p.MAX_A+p.GRID_A_SHIFT),p.N_A_NS)';
grid_a_NS = exp(grid_a_preshift_NS) - p.GRID_A_SHIFT;

% Discretize the exogenous stochastic process
[grid_z prob_z] = Rouwenhorst(p.RHO_Z,p.VAR_Z_INN,p.N_Z);
grid_z = exp(grid_z);
[v d] = eig(prob_z');
s_z = v(:,1)/sum(v(:,1));

% Set interest rates and wages to the steady-state levels.
RToday = R_SS;
wToday = W_SS;

% bisection
bisec_min = 0.75;
bisec_max = 0.995;

meanY  = 1.0; % The initial guess of Yss.
beta = (bisec_min + bisec_max)/2;

% Initial guess for the policy functions for consumption and labor supply

pf_c = zeros(p.N_A,p.N_Z);
pf_n = zeros(p.N_A,p.N_Z);
pf_sav = zeros(p.N_A,p.N_Z);
values = zeros(p.N_A,p.N_Z);

for i_a = 1:p.N_A

    aToday = grid_a(i_a);

    for i_z = 1:p.N_Z
    
        zToday = grid_z(i_z);
        pf_c(i_a,i_z) = 0.3 + 0.1*aToday;
        pf_n(i_a,i_z) = (wToday*zToday*pf_c(i_a,i_z)^(-p.GAMMA))^(1.0/p.PSI);
        
        % Initial guess for the value function is for interpolating it below....
        if (p.GAMMA == 1.0)
            values(i_a,i_z) = log(aToday + zToday)/(1.0-beta);
        else
            values(i_a,i_z) = ((aToday + zToday)^(1.0-p.GAMMA)/(1.0-p.GAMMA))/(1.0-beta);
        end

    end

end

t=cputime;

for beta_iter = 1:p.MAX_BETA_ITER

    beta = (bisec_min + bisec_max)/2

    tauToday = 1.0/s_z(3)*p.B4YSS*4.0*meanY*(RToday-1.0)/RToday/p.TAU_BAR(3);
    dToday   = meanY*(1.0 - wToday);

    [pf_c pf_n pf_sav] = HH_opt_EGM(beta,p,grid_a,grid_z,prob_z,wToday,RToday,tauToday,dToday,pf_c,pf_n,pf_sav);

    [dist meanC meanN meanA] = HH_dist(beta,p,grid_a,grid_a_NS,grid_z,prob_z,wToday,RToday,tauToday,dToday,pf_c,pf_n,pf_sav);

    disp([beta_iter meanC meanN meanA]);
    meanY = meanC;

    % Evaluate convergence
    if ( abs(meanA/(4.0*meanY) - p.B4YSS ) < p.BETA_ERR_TOL ); break; end;

    % Update the guess using the bisection method
    if ( meanA/(4.0*meanY) > p.B4YSS )
        bisec_max = beta;
    else
        bisec_min = beta;
    end

end

disp(sprintf('Total elapsed time = %6.5f seconds',cputime-t));

save HANK_SS grid_a_NS pf_c pf_n pf_sav dist;

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

