%% Example: price setting without feedback
% This example solves the simple model in Afrouzi and Yang (2020). 

clear;
addpath('../src');
%% parameters 

% DRIP parameters
omega   = sqrt(0.000045);
beta    = 0.99;
su      = 0.01;
phiL    = 1;
phiH    = 1.5;

% other
L       = 20; % length of irfs

%% Solvd steady state 
solL       = Drip(omega,beta,1,su/phiL,1);
solH       = Drip(omega,beta,1,su/phiH,1);
%% Aftermath of more hawkish policy
solH.tr   = Trip(solH,solL.ss.Sigma_1,'T',40);

% get irfs under the steady state information structure 
solH.irfs.ss = irfs(solH,'length',L);

% get irfs under the transition path information structure 
solH.irfs.tr = irfs(solH,'length',L,'trip',solH.tr);

%% Aftermath of more dovish policy 
solL.tr   = Trip(solL,solH.ss.Sigma_1,'T',40);

% get irfs under the steady state information structure 
solL.irfs.ss = irfs(solL,'length',L);

% get irfs under the transition path information structure 
solL.irfs.tr = irfs(solL,'length',L,'trip',solL.tr);
%% Plot irfs of fundamental and price 
figure;
subplot(2,2,[1,2]);
hold on;
box on;
grid on;
xlim([1,ex1.irfs.ss.T]);
title('IRFs to 1 Std. Dev. Expansionary Shock','interpreter','latex', ...
    'fontsize',12);
plot(1:ex1.irfs.ss.T,reshape(ex1.irfs.ss.x(1,1,:),[ex1.irfs.ss.T,1]), ...
    'LineWidth',3,'Color','k')
plot(1:ex1.irfs.ss.T,reshape(ex1.irfs.ss.a(1,1,:),[ex1.irfs.ss.T,1]), ...
    'LineWidth',3,'Color','b')
plot(1:ex1.irfs.tr.T,reshape(ex1.irfs.tr.a(1,1,:),[ex1.irfs.tr.T,1]), ...
    'LineWidth',3,'Color','r')
legend({'Nominal Agg. Demand ($q$)','Price IRF in Steady State ($p$)', ...
    'Price IRF in Transition Path ($p$)'}, ...
    'interpreter','latex','fontsize',10,'location','southeast')

%% Plot irfs of output 
% define output as the difference between nominal GDP and price 
y_ss = reshape(ex1.irfs.ss.x(1,1,:)-ex1.irfs.ss.a(1,1,:),[ex1.irfs.ss.T,1]);
y_tr = reshape(ex1.irfs.tr.x(1,1,:)-ex1.irfs.tr.a(1,1,:),[ex1.irfs.tr.T,1]);

subplot(2,2,3);
hold on;
box on;
grid on;

xlim([1,ex1.irfs.ss.T]);
title('Output','interpreter','latex', ...
    'fontsize',14);
plot(1:ex1.irfs.ss.T,y_ss, ...
    'LineWidth',3,'Color','b')
plot(1:ex1.irfs.ss.T,y_tr, ...
    'LineWidth',3,'Color','r')
legend({'Steady State','Transition Path'}, ...
    'interpreter','latex','fontsize',10,'location','northeast')

%% Plot irfs of inflation 
% define inflation as the growth in price
pi_ss = reshape(ex1.irfs.ss.a(1,1,2:end)-ex1.irfs.ss.a(1,1,1:end-1),[ex1.irfs.ss.T-1,1]);
pi_tr = reshape(ex1.irfs.tr.a(1,1,2:end)-ex1.irfs.tr.a(1,1,1:end-1),[ex1.irfs.tr.T-1,1]);
pi_ss = [ex1.irfs.ss.a(1,1,1);pi_ss];
pi_tr = [ex1.irfs.tr.a(1,1,1);pi_tr];

subplot(2,2,4);
hold on;
box on;
grid on;
xlim([1,ex1.irfs.ss.T]);
title('Inflation','interpreter','latex', ...
    'fontsize',14);
plot(1:ex1.irfs.ss.T,pi_ss, ...
    'LineWidth',3,'Color','b')
plot(1:ex1.irfs.ss.T,pi_tr, ...
    'LineWidth',3,'Color','r')
legend({'Steady State','Transition Path'}, ...
    'interpreter','latex','fontsize',10,'location','northeast')
