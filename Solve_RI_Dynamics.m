%% Afrouzi and Yang (2019)

%% Inputs
% omega     : Cost of information
% beta	    : Discount factor
% A 	    : x=Ax+Qu
% Q 	    : x=Ax+Qu
% H	    : v=-0.5(a'-x'H)(a-H'x)
% Omega_init: Initial guess for benefit matrix Omega
% Sigma_init: Initial guess for prior covariance matrix 

%% Outputs
% Sigma_1   : Steady-state Prior uncertainty
% Sigma_p   : Steady-state posterior uncertainty
% Lambda    : Shadow matrix on the no-forgetting constriant -- Lambda*(Sigma-Sigma_p) = 0 
% Omega     : Dynamic benefit matrix 
% Y         : Weight vector for evolution of actions
% Sigma_z   : Covariance matrix of the rational inattention error
% K	    : Kalman gain matrix



function [Sigma1,Sigma_p,Lambda,Omega0,Y,Sigma_z,K] = ...
                    Solve_RI_Dynamics(omega,beta,A,Q,H,Omega_init,Sigma_init) 

    tol_err = 1e-8      ;   % tolerance level
    w       = 0.1       ;   % update weight

    Omega_c = H*H'      ;
    
    [n,m]   = size(H)   ;
    I       = eye(n,n)  ;

    % Guess Sigma_t|t-1 and Omega
    Omega0  = Omega_init ;
    Sigma0  = Sigma_init ;

    err     = 1 ;
    iter    = 0 ;
    
	SqRSigma    = sqrtm(Sigma0);
    
    % Loop
    while err > tol_err
		
		S_Om_S  = SqRSigma*Omega0*SqRSigma ;
        [U,D] = eig(S_Om_S) ;
        D = diag((abs(diag(D))>1e-10).*diag(D)) ;
        U = (abs(U)>1e-10).*U ;

        
        Lambda      = U*max(omega*I - D,0)*U'; 
        Sigma_p     = omega*SqRSigma*U/(max(D,omega*I))*U'*SqRSigma;
    
        Sigma1      = A*Sigma_p*A' + Q*Q' ;
        err         = norm(Sigma1 - Sigma0,'fro')/sqrt(n);

        Sigma0      = w*Sigma1 + (1-w)*Sigma0 ;
        
        SqRSigma    = sqrtm(Sigma0);
        invSqRSigma = pinv(SqRSigma);
        
        Omega0      = Omega_c + beta*A'*invSqRSigma*(omega*I - Lambda) ... 
					  *invSqRSigma*A   ;

        SqRSigma = (abs(SqRSigma)>1e-10).*SqRSigma ;
        Omega0 = (abs(Omega0)>1e-10).*Omega0 ;

        iter        = iter + 1;
    end
   
    inv_Sigma1  = pinv(Sigma1) ;
    Y = (I - Sigma_p*inv_Sigma1)'*H ;
    Sigma_z = H'*(Sigma_p - Sigma_p*inv_Sigma1*Sigma_p)*H ;
    K = Sigma1*Y*pinv(Y'*Sigma1*Y + Sigma_z) ;
    
end
