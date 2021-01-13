function [y_signal,Sigma,Omega_bar,D_eid,K] = Solve_RI_dynamics_SignalBase(beta,omega,w,A,Q,SigmaInit,OmegaInit)
%This function solves dynamic rational inattention problems. 
%beta is the discount factor of the agent.
%omega is the marginal cost of information
%the fundamental follows a process w'*U(t) where U(t) is the state of the  
%economy with the state space respresentation U(t)=A*U(t-1)+Q*u(t), Q=b*b'. 
%u(t) is the innovation to the state. 

%% Parameters

T       =length(w)  ;   %length of the state space
I       =eye(T)     ;   %definig the identity matrix
Omega   =w*w'       ;   
slow_inc=0.8        ;   %weight on old guess when updating guesses
iter    =0          ;
tol_y   =1          ;

% Initial Guess
Omega_bar   = OmegaInit ; %initial guess for Omega
Sigma0      = SigmaInit ; %initial guess for Sigma

while tol_y>1e-6
    
    SqRSigma    = sqrtm(Sigma0)         ;
    invSqRSigma = pinv(SqRSigma)        ;
    S_Om_S      = SqRSigma*Omega_bar*SqRSigma ;
    [U,D]       = eig(S_Om_S) ;
    
    d   = diag(D)     ;
    
    ii = find(d>omega) ;
    y   = zeros(T,length(ii)) ;
	
    temp1= zeros(T,T) ;
    temp2= zeros(T,T) ;
    for i = 1:length(ii)
        y(:,i) =invSqRSigma*U(:,ii(i)) ;
        temp1 = temp1 + (d(ii(i))-omega)*y(:,i)*y(:,i)' ; 
        temp2 = temp2 + (1-omega/d(ii(i)))*y(:,i)*y(:,i)' ; 
    end
	
	% Update
    Omega_bar   = Omega - beta*A'*(temp1)*A + beta*A'*Omega_bar*A ;
    Sigma1      = A*(Sigma0 - Sigma0*(temp2)*Sigma0)*A' + Q*Q' ;
    
    tol_y       = norm(Sigma1 - Sigma0,'fro')/sqrt(T) ;    
	
    Sigma0      = slow_inc*Sigma1 + (1-slow_inc)*Sigma0;  

    iter = iter + 1 ;
    if mod(iter,5000) == 0
        fprintf('      iter = %d,  errA = %10.7f \n', iter, tol_y)    
    end
    
    if iter > 100000 
        break ;
    end
end
%fprintf('      iteration done = %d,  i(d>omega) = %d,  Evalue = %5.2g \n', iter, ii, d(ii))    

%% retrun outputs
Sigma   =Sigma0;
K       =zeros(T,T) ;
y_signal=zeros(T,length(ii)) ;
for i = 1:length(ii)
    K = K + (1-omega/d(ii(i)))*Sigma*y(:,i)*y(:,i)' ; 
    y_signal(:,i)=sqrt(1-omega/d(ii(i)))*y(:,i) ;
end
D_eid = d ;

end

