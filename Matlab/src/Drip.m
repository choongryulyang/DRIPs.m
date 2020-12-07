%% Afrouzi and Yang (2019)
%% This code solves for the steady state of the dynamic 
%% multivariate rational inattention problem in LQG setting.
%
%     Drip(omega,beta,A,Q,H, kwargs...) -> Drip
%
% Solves for the steady state of a Dynamic Rational Inattention Problem defined 
% by the arguments and stores the solution in a Drip type.
% See [Afrouzi and  Yang (2020)](http://afrouzi.com/dynamic_inattention.pdf) for details.
%
%% ARGUMENSTS
%
% The function takes the primitives of the Drip as arguments:
%     * omega  : Cost of information
%     * beta   : Discount factor
%     * A      : Transition matrix: x=Ax+Qu
%     * Q      : Std. Dev. matrix: x=Ax+Qu
%     * H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)
%
%% OPTIONAL ARGUMENSTS (kwargs...)
%
% Default values are set unless specified otherwise by user.
%     * fcap      [= false]    [if `true` then solves the problem with fixed capacity = ω bits]
%     * initOmega [= H*H']     [initial guess for steady state information matrix]
%     * initSigma [= A*A'+Q*Q'][initial guess for steady state prior]
%     * w         [= 1]        [updating weight on the new guess in iteration]
%     * tol       [= 1e-4]     [tolerance level for convergence]
%     * maxit     [= 10000]    [maximum number of iterations]
%
%% OUTPUTS
%  The function returns a structure with the primitives and the solution objects:
%     * ss : a structure with the solution to the steady state of the DRIP with the following fields:
%           * Y       : Weight vector for evolution of actions
%           * Sigma_z : Covariance matrix of the rational inattention error
%           * K       : Kalman gain matrix
%           * Sigma_1 : Steady-state prior covariance matrix under the solution
%           * Sigma_p : Steady-state posterior covariance matrix under the solution
%           * Omega   : Dynamic benefit matrix
%     * also contains fields for the primitives of the problem (omega,beta,A,Q,H)
%
%% EXAMPLE 
% >> p = Drip(omega,beta,A,Q,H)
% 

function p = Drip(omega,beta,A,Q,H,varargin) 

    % parse optional inputs
    args = inputParser;
    addOptional(args,'initOmega',H*H');
    addOptional(args,'initSigma',A*A'+Q*Q');
    addOptional(args,'tol',1e-4);
    addOptional(args,'fcap',false,@islogical);
    addOptional(args,'w',1);
    addOptional(args,'maxit',10000);

    parse(args,varargin{:});

    % initialize 

    [n,~]      = size(H)   ;                % get the dimension of the state n 
    p.ss.err   = 1 ;                        % initialize convergence error
    iter       = 0 ;                        % initialize number of iterations
    Sigma0     = args.Results.initSigma ;   % initial guess for Omega
    p.ss.Omega = args.Results.initOmega ;   % initial guess for Sigma
    
    % other variables to be used in the code:
    SqRSigma   = sqrtm(Sigma0);
    Omega_c    = H*H'      ;
    I          = eye(n)    ;

    % if the user wants solution with fixed capacity, set kappa = omega
    if args.Results.fcap == true 
        kappa     = omega ;
    end 	
    
    % Loop (stop if converged or iteration reached maxit)
    while (p.ss.err > args.Results.tol) && (iter <= args.Results.maxit)
		
        % get the eigenvalue and eigen vectors of SqSigma*Omega*SqSigma
		X     = SqRSigma*p.ss.Omega*SqRSigma ;
        [U,D] = eig(X) ;
        D = real(D);
        U = real(U);
        
        % if the user wants solution with fixed capacity, find the
        % corresponding lagrange multiplier on the capacity constriant:
        if args.Results.fcap == true
            omega = (2^(2*kappa)/det(max(D,omega*I)))^(-1/n);
        end
        
        % use the policy function in Afrouzi and Yang to get the implied
        % posterior covariance:
        p.ss.Sigma_p = omega*SqRSigma*U/(max(D,omega*I))*U'*SqRSigma;
    
        % get the new guess for prior covariance and calculate the
        % convergence error:
        p.ss.Sigma_1 = A*p.ss.Sigma_p*A' + Q*Q' ;
        p.ss.err     = norm(p.ss.Sigma_1 - Sigma0)/norm(Sigma0);

        % update the guess for prior covariance given update weight w
        Sigma0      = args.Results.w*p.ss.Sigma_1 + (1-args.Results.w)*Sigma0 ;
        
        % use the Euler equation in Afrouzi and Yang to update guess for
        % Omega given update weight w:
        SqRSigma    = real(sqrtm(Sigma0));
        invSqRSigma = real(pinv(SqRSigma));
        p.ss.Omega   = args.Results.w ...
                        *(Omega_c + beta*A'*invSqRSigma*U*(min(D,omega*I))*U'*invSqRSigma*A) ...
                        +(1-args.Results.w)*p.ss.Omega;
        p.ss.Omega = (abs(p.ss.Omega)>1e-10).*p.ss.Omega ;

        % update number of iterations
        iter        = iter + 1;
    end
    
    % produce matrices Y, K, Sigma_z as in Afrouzi and Yang for output
    inv_Sigma1  = pinv(p.ss.Sigma_1) ;
    p.ss.Y       = (I - p.ss.Sigma_p*inv_Sigma1)'*H ;
    p.ss.Sigma_z = H'*(p.ss.Sigma_p - p.ss.Sigma_p*inv_Sigma1*p.ss.Sigma_p)*H ;
    p.ss.K       = p.ss.Sigma_1*p.ss.Y*pinv(p.ss.Y'*p.ss.Sigma_1*p.ss.Y + p.ss.Sigma_z) ;
    
    % store the primitives of the problem
    p.omega = omega;
    p.beta = beta;
    p.A = A;
    p.Q = Q;
    p.H = H;
end
