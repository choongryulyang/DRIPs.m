%% Afrouzi and Yang (2019)
%% This code solves for the steady state of the dynamic 
%% multivariate rational inattention problem in LQG setting.

%     Drip(ω,β,A,Q,H, kwargs...) -> Drip
%
% Solves for the steady state of a Dynamic Rational Inattention Problem defined 
% by the arguments and stores the solution in a Drip type.
% See [Afrouzi and  Yang (2020)](http://afrouzi.com/dynamic_inattention.pdf) for details.
%
%% ARGUMENSTS
%
% The function takes the primitives of the Drip as arguments:
%     * ω      : Cost of information
%     * β      : Discount factor
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
%     * Y       : Weight vector for evolution of actions
%     * Sigma_z : Covariance matrix of the rational inattention error
%     * K       : Kalman gain matrix
%     * Sigma_1 : Steady-state prior covariance matrix under the solution
%     * Sigma_p : Steady-state posterior covariance matrix under the solution
%     * Omega   : Dynamic benefit matrix
%
%% EXAMPLE 
% >> P = Drip(ω,β,A,Q,H)
% 




function sol = Drip(omega,beta,A,Q,H,varargin) 

    % parse optional inputs 
    
    p = inputParser;
    addOptional(p,'initOmega',H*H');
    addOptional(p,'initSigma',A*A'+Q*Q');
    addOptional(p,'tol',1e-4);
    addOptional(p,'fcap',false,@islogical);
    addOptional(p,'w',1);
    addOptional(p,'maxit',10000);

    parse(p,varargin{:});

    % initialize 

    [n,~]     = size(H)   ;
    I         = eye(n)    ;
    sol.err   = 1 ;
    iter      = 0 ;
    Sigma0    = p.Results.initSigma ;
    sol.Omega = p.Results.initOmega ;
    SqRSigma  = sqrtm(Sigma0);
    Omega_c   = H*H'      ;

    if p.Results.fcap == true 
        kappa     = omega ;
    end 	
    
    % Loop
    while sol.err > p.Results.tol
		
		S_Om_S  = SqRSigma*sol.Omega*SqRSigma ;
        [U,D] = eig(S_Om_S) ;
        D = real(D);
        U = real(U);

        if p.Results.fcap == true
            omega = (2^(2*kappa)/det(max(D,omega*I)))^(-1/n);
        end
        
        sol.Sigma_p = omega*SqRSigma*U/(max(D,omega*I))*U'*SqRSigma;
    
        sol.Sigma_1 = A*sol.Sigma_p*A' + Q*Q' ;
        sol.err     = norm(sol.Sigma_1 - Sigma0)/norm(Sigma0);

        Sigma0      = p.Results.w*sol.Sigma_1 + (1-p.Results.w)*Sigma0 ;
        
        SqRSigma    = real(sqrtm(Sigma0));
        invSqRSigma = real(pinv(SqRSigma));
        
        sol.Omega   = p.Results.w ...
                        *(Omega_c + beta*A'*invSqRSigma*U*(min(D,omega*I))*U'*invSqRSigma*A) ...
                        +(1-p.Results.w)*sol.Omega;

        sol.Omega = (abs(sol.Omega)>1e-10).*sol.Omega ;

        iter        = iter + 1;
    end
   
    inv_Sigma1  = pinv(sol.Sigma_1) ;
    sol.Y       = (I - sol.Sigma_p*inv_Sigma1)'*H ;
    sol.Sigma_z = H'*(sol.Sigma_p - sol.Sigma_p*inv_Sigma1*sol.Sigma_p)*H ;
    sol.K       = sol.Sigma_1*sol.Y*pinv(sol.Y'*sol.Sigma_1*sol.Y + sol.Sigma_z) ;
    
end
