%          Trip(p::Drip,Sigma0, kwargs...) -> Trip
%
%     Solves for the transition dynamics of the optimal information structure starting
%     from the initial prior distribution with covariance matrix `Sigma0`.
%     See [Afrouzi and  Yang (2019)](http://afrouzi.com/dynamic_inattention.pdf)
%     for details.
%
%% ARGUMENTS
%     * p       : solution to the steady state of Drip (output of Drip.m function)
%     * Sigma0  : initial prior covariance matrix
%
%% OPTIONAL ARGUMENTS
%     * T     = 100  : guess for convergence time until steady state
%     * tol   = 1e-4 : tolerance for convergence
%     * maxit = 1000 : max iterations
%
%% OUTPUTS
%     Returns a structure with the following fields for the transition path of
%     the optimal information structure:
%
%     Sigma_1s : sequence of prior covariance matrices
%     Sigma_ps : sequence of posterior covariance matrices
%     Omegas   : sequence of information benefit matrices
%     Ds       : eigenvalues of Sigma_t^(0.5)Omega_tSigma_t^(0.5) over time (marginal values of information)
%     err      : convergence err
%     con_err  : distance of Sigma_T from steady state prior
% 
%% EXAMPLE
% >> p = solve_drip(ω,β,A,Q,H)
% >> Sigma0 = 0.1*p.Sigma_1;
% >> pt = Trip(p,Sigma0);

function sol = Trip(p,Sigma0,varargin) 

    % parse optional inputs 
    args = inputParser;
    addOptional(args,'T',100);
    addOptional(args,'tol',1e-4);
    addOptional(args,'maxit',1000);

    parse(args,varargin{:});

    % initialize 
    Omegas0     = repmat(p.ss.Omega,1,1,args.Results.T);
    sol.Omegas  = Omegas0;

    Sigma_1s0    = repmat(p.ss.Sigma_1,1,1,args.Results.T); Sigma_1s0(:,:,1) = Sigma0;
    sol.Sigma_1s = Sigma_1s0;

    [n,~]     = size(p.H)   ;

    sol.Ds    = zeros(n,args.Results.T);

    sol.Sigma_ps = repmat(p.ss.Sigma_p,1,1,args.Results.T);

    iter = 0;
    sol.err  = 1;
    I    = eye(n);

    % loop
    while (sol.err > args.Results.tol) && (iter <= args.Results.maxit)
       % Given Omegas, find Sigmas using the law of motion for priors
      for i = 1:1:args.Results.T-1
        SqSigma = real(sqrtm(Sigma_1s0(:,:,i)));
        [U, D]  = eig(SqSigma*Omegas0(:,:,i)*SqSigma);
        U       = real(U);
        D       = real(D);
        sol.Ds(:,i)           = sort(diag(D));
        sol.Sigma_ps(:,:,i)   = p.omega*SqSigma*U/(max(D,p.omega*I))*U'*SqSigma;
        sol.Sigma_1s(:,:,i+1) = p.Q*p.Q' + p.A*sol.Sigma_ps(:,:,i)*p.A';
      end
      % Given Sigma_1s, find Omegas using the Euler equation
      for i = args.Results.T-1:-1:1
          SqSigma      = real(sqrtm(sol.Sigma_1s(:,:,i+1)));
          invSqSigma   = real(pinv(SqSigma));
          [U, D]       = eig(SqSigma*Omegas0(:,:,i+1)*SqSigma);
          D            = real(D);
          U            = real(U);
          sol.Omegas(:,:,i)   = p.H*p.H' + p.beta*p.A'*invSqSigma*U* ...
                                min(D,p.omega*I)*U'*invSqSigma*p.A;
      end
      sol.err = norm(sol.Sigma_1s(:)-Sigma_1s0(:))/norm(Sigma_1s0(:));
      Sigma_1s0 = real(sol.Sigma_1s);
      Omegas0 = real(sol.Omegas);
      iter = iter + 1;
    end
    sol.con_err = norm(p.Q*p.Q' + p.A*sol.Sigma_ps(:,:,end-1)*p.A'-p.ss.Sigma_1)/norm(p.ss.Sigma_1);
    SqSigma       = real(sqrtm(sol.Sigma_1s(:,:,end)));
    [~, D]        = eig(SqSigma*sol.Omegas(:,:,end)*SqSigma);
    sol.Ds(:,end) = sort(diag(real(D)));

  end 