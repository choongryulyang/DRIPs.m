
% Sigma_1   : Steady-state Prior uncertainty
% Sigma_p   : Steady-state posterior uncertainty
% Lambda    : Lambda*(Sigma-Sigma_p) = 0
% Omega     : Forward looking component of the FOC
% Y         : Weight vector for evolution of actions
% Sigma_z   : Covariance matrix of the rational inattention error

function [Sigma1,Sigma_p,Lambda,Omega0,Y,Sigma_z,K,X] = ...
                    Solve_RI_Dynamics(phi,beta,A,Q,H,Omega_init,Sigma_init) 

    tol_err = 1e-8     ;   % tolerance level
    w       = 0.5       ;   % update weight

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

        %D = diag((abs(diag(D))>1e-8).*diag(D)) ;
        
        Lambda      = U*max(phi*I - D,0)*U'; 
        Sigma_p     = phi*SqRSigma*U*inv(max(D,phi*I))*U'*SqRSigma;
    
        Sigma1      = A*Sigma_p*A' + Q*Q' ;
        err         = norm(Sigma1 - Sigma0,'fro')/n;

        Sigma0      = w*Sigma1 + (1-w)*Sigma0 ;
        
        SqRSigma    = sqrtm(Sigma0);
        invSqRSigma = pinv(SqRSigma);
        
        Omega0      = Omega_c + beta*A'*invSqRSigma*(phi*I - Lambda) ... 
					  *invSqRSigma*A   ;
        
        iter        = iter + 1;
    end

    Y = (I - Sigma_p*pinv(Sigma1))'*H ;
    Sigma_z = H'*(Sigma_p - Sigma_p*pinv(Sigma1)*Sigma_p)*H ;
    X           = Sigma1 - Sigma_p;
    % Add Kalman gain matrix
    % Notice I use pseudo-inverse matrix since sometimes Sigma_z is
    % singluar
    K = Sigma1*Y*pinv(Y'*Sigma1*Y + Sigma_z) ;
    
end