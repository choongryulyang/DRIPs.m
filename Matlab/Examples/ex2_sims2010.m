%% Example: Sims 2010
% This code solves for the Sims 2010 and its extesion in Afrouzi and Yang
% (2020)
clear;
addpath('../src');
%% parameters 
    % DRIP parameters
    omega = 1;
    beta  = 0.9;
    A     = [0.95 0; 0 0.4];
    Q     = sqrt([0.0975 0; 0 0.84]);
    H     = [1; 1];

%% Solution and Performance for Benchmark Parameterization
    % Solve and display the optimal posterior covariance matrix:
    bp    = Drip(omega,beta,A,Q,H);
    disp('Posterior covariance matrix in the steady state:')
    disp(bp.ss.Sigma_p)

    % Measure performance
    disp('Median time for obtaining the solution (in seconds):')
    disp(timeit(@() Drip(omega,beta,A,Q,H)))

%% Solution and Performance for Transition Dynamics
    % Define initial prior as 0.01*Steady state prior
    Sigma0 = 0.01*bp.ss.Sigma_1;
    
    % Solve for the transition dynamics and store in bp.tr
    bp.tr = Trip(bp,Sigma0);

    % Measure performance
    disp('Median time for obtaining the solution for transition dynamics (in seconds):')
    disp(timeit(@() Trip(bp,Sigma0,'T',20))); % 'T' option specifies time until convergence 

%% Plot marginal value of information
    figure;
    hold on;
    box on;
    grid on;
    xlim([1,15])
    title('Marginal Value of Information on the Transition Path', ...
        'interpreter','latex','fontsize',12);
    plot(1:15,bp.tr.Ds(1,1:15), ...
        'LineWidth',3,'Color','b')
    plot(1:15,bp.tr.Ds(2,1:15), ...
        'LineWidth',3,'Color','r')
    plot(1:15,bp.omega*ones(1,15),'--k','LineWidth',2)
    legend({'High eigenvalue ($d_1$)','Low eigenvalue ($d_2$)', ...
        'Marginal cost ($\omega$)'}, ...
        'interpreter','latex','fontsize',10,'location','east')
    
    
%% Calculate irfs 
    % Get irfs for steady state info structure and store in bp.irfs.ss
    bp.irfs.ss = irfs(bp,'length',20);
    
    % Get irfs for transition dynamics and store in bp.irfs.tr
    bp.irfs.tr = irfs(bp,'trip',bp.tr,'length',20);

%% Plot IRFs
    figure;
    
    % irfs to slow moving shock 
    subplot(2,1,1);
    hold on;
    box on;
    grid on;
    xlim([1,bp.irfs.ss.T]);
    title('IRFs to Slow Moving Shock','interpreter','latex', ...
        'fontsize',12);
    plot(1:bp.irfs.ss.T,reshape(bp.irfs.ss.x(1,1,:),[bp.irfs.ss.T,1]), ...
        'LineWidth',3,'Color','k')
    plot(1:bp.irfs.ss.T,reshape(bp.irfs.ss.a(1,1,:),[bp.irfs.ss.T,1]), ...
        'LineWidth',3,'Color','b')
    plot(1:bp.irfs.tr.T,reshape(bp.irfs.tr.a(1,1,:),[bp.irfs.tr.T,1]), ...
         'LineWidth',3,'Color','r')
    legend({'Shock','Price IRF in Steady State ($p$)', ...
        'Price IRF in Transition Path ($p$)'}, ...
        'interpreter','latex','fontsize',10,'location','northeast')

    % irfs to fast moving shock 
    subplot(2,1,2);
    hold on;
    box on;
    grid on;
    xlim([1,bp.irfs.ss.T]);
    title('IRFs to Slow Moving Shock','interpreter','latex', ...
        'fontsize',12);
    plot(1:bp.irfs.ss.T,reshape(bp.irfs.ss.x(2,2,:),[bp.irfs.ss.T,1]), ...
        'LineWidth',3,'Color','k')
    plot(1:bp.irfs.ss.T,reshape(bp.irfs.ss.a(1,2,:),[bp.irfs.ss.T,1]), ...
        'LineWidth',3,'Color','b')
    plot(1:bp.irfs.tr.T,reshape(bp.irfs.tr.a(1,2,:),[bp.irfs.tr.T,1]), ...
         'LineWidth',3,'Color','r')
    legend({'Shock','Price IRF in Steady State ($p$)', ...
        'Price IRF in Transition Path ($p$)'}, ...
        'interpreter','latex','fontsize',10,'location','northeast')
