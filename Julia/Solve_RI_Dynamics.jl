## Afrouzi and Yang (2019)
## This code solves for the steady state of the dynamic 
## multivariate rational inattention problem in LQG setting.
##
## Julia code written by Miguel Acosta / Modified by Hassan Afrouzi

## Inputs
# ω         : Cost of information
# β     : Discount factor
# A         : x=Ax+Qu
# Q         : x=Ax+Qu
# H     : v=-0.5(a'-x'H)(a-H'x)
# Ωguess    : Initial guess for benefit matrix Ω
# Σguess    : Initial guess for prior covariance matrix 

## Outputs
# Σ_1       : Steady-state Prior uncertainty
# Σ_p       : Steady-state posterior uncertainty
# Λ         : Shadow matrix on the no-forgetting constriant -- Λ*(Σ-Σ_p) = 0 
# Ω         : Dynamic benefit matrix 
# Y         : Weight vector for evolution of actions
# Σ_z       : Covariance matrix of the rational inattention error
# K     : Kalman gain matrix

using LinearAlgebra

struct drip
    ω; β; A; Q; H;          #Primitives
    K; Y; Σ_z; Σ_p; Σ_1; Ω; #Solution
end

struct drip_init
    Σ::Array{Float64}
    Ω::Array{Float64}
end

struct dripirfs
    T::Int
    x::Array{Float64}
    x_hat::Array{Float64}
    a::Array{Float64}
end

 
function Solve_RI_Dynamics(ω,β,A,Q,H;Ω0=H*H',Σ0=A*A'+Q*Q')
    tol_err = 1e-8       # tolerance level
    w       = 0.1        # update weight
    maxit   = 10000

    Ω_c = H*H'
    (n,m)     = length(size(H)) == 2 ? size(H) : (size(H,1),1)
    eye       = Matrix{Float64}(I,n,n)

    err     = 1
    iter    = 0
    
    SqRΣ    = sqrt(Σ0)

    Σ1  = Matrix{Float64}(I,n,n)
    Σ_p = Matrix{Float64}(I,n,n)
    Λ   = Matrix{Float64}(I,n,n)
    
    # Loop
    while (err > tol_err) & (iter < maxit)
        
    S_Om_S  = SqRΣ*Ω0*SqRΣ
        U       = eigvecs(S_Om_S)
        D       = eigvals(S_Om_S)
        D = getreal(D)
        D = diagm(abs.(D).>1e-10).*diagm(D) + (diagm(abs.(D).<=1e-10))*1e-8
        U = getreal(U) #(abs.(U).>1e-10).*U ;
                  
        Λ       = U*max.(ω*eye - D,0.0)*U'; 
        Σ_p     = ω*SqRΣ*U/(max.(D,ω*eye))*U'*SqRΣ;
        Σ_p     = getreal(Σ_p)        

        
        Σ1      = A*Σ_p*A' + Q*Q'
        err         = norm(Σ1 - Σ0,2)/sqrt(n)
        
        Σ0      = w*Σ1 + (1-w)*Σ0
        
        SqRΣ    = sqrt(Σ0);
        invSqRΣ = pinv(SqRΣ);
        
        Ω0      = w*(Ω_c + β*A'*invSqRΣ*(ω*eye - Λ)*invSqRΣ*A)+(1-w)*Ω0;

        SqRΣ = (abs.(SqRΣ).>1e-10).*SqRΣ + diagm(abs.(diag(SqRΣ)).<=1e-10)*1e-8
        Ω0   = (abs.(Ω0).>1e-10).*Ω0     + diagm(abs.(diag(Ω0)).<=1e-10)*1e-8

        iter += 1
    end

    if iter == maxit
        print("RI Code hit maxit\n")
    end

    Σ1      = (abs.(Σ1).>1e-10).*Σ1 + diagm(abs.(diag(Σ1)).<=1e-10)*1e-8
    inv_Σ1  = pinv(Σ1) ;
    Y = (eye - Σ_p*inv_Σ1)'*H ; 
    Σ_z = H'*(Σ_p - Σ_p*inv_Σ1*Σ_p)*H ;
    K = Σ1*Y*pinv(Y'*Σ1*Y + Σ_z) ;

    sol=drip(ω,β,A,Q,H,K,Y,Σ_z,Σ_p,Σ1,Ω0)
    return(sol)
end



########## Aux. Functions ############
function getreal(M)
    if maximum(abs.(imag.(M))) < 1e-10
        return(real.(M))
    else
        print("Your matrix has complex stuff")
    end
end


function infinitesum(func; tol = 1e-8,maxit = 1000,start=0)
    diff  = 1.0
    infsum = func(start)
    it    = start + 1
    while (diff > tol) & (it < maxit)
        func_it = func(it)
        infsum += func_it
        diff = maximum(func_it)
        it += 1
    end
    return(infsum)
end

function dripirfs(P::drip,T::Int)
    (n,m) = length(size(P.H)) == 2 ? size(P.H) : (size(P.H,1),1)
    (_,k) = length(size(P.Q)) == 2 ? size(P.Q) : (size(P.Q,1),1)
    x     = zeros(n,k,T);
    x_hat = zeros(n,k,T);
    a     = zeros(m,k,T);
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1 
                x[:,kk,ii]     = P.Q*e_k;
                x_hat[:,kk,ii] = (P.K*P.Y')*(x[:,kk,ii]);
            else 
                x[:,kk,ii]     =P.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =P.A*x_hat[:,kk,ii-1]+(P.K*P.Y')*(x[:,kk,ii]-P.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= P.H'*x_hat[:,kk,ii];
        end
    end
    return(dripirfs(T,x,x_hat,a))
end


