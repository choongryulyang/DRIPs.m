## Afrouzi and Yang (2019)
## This code solves for the steady state of the dynamic 
## multivariate rational inattention problem in LQG setting.
##
## Julia code written by Miguel Acosta

## Inputs
# ω         : Cost of information
# β	    : Discount factor
# A 	    : x=Ax+Qu
# Q 	    : x=Ax+Qu
# H	    : v=-0.5(a'-x'H)(a-H'x)
# Ωguess    : Initial guess for benefit matrix Ω
# Σguess    : Initial guess for prior covariance matrix 

## Outputs
# Σ_1       : Steady-state Prior uncertainty
# Σ_p       : Steady-state posterior uncertainty
# Λ         : Shadow matrix on the no-forgetting constriant -- Λ*(Σ-Σ_p) = 0 
# Ω         : Dynamic benefit matrix 
# Y         : Weight vector for evolution of actions
# Σ_z       : Covariance matrix of the rational inattention error
# K	    : Kalman gain matrix


using LinearAlgebra

function Solve_RI_Dynamics(ω::Float64,β::Float64,
                           A::Array{Float64},Q::Array{Float64},H::Array{Float64},
                           Ωguess::Array{Float64},Σguess::Array{Float64})

    tol_err = 1e-8       # tolerance level
    w       = 0.1        # update weight
    maxit   = 1e6

    Ω_c = H*H'
    (n,m)     = length(size(H)) == 2 ? size(H) : (size(H,1),1)
    eye       = Matrix{Float64}(I,n,n)

    # Guess Σ_t|t-1 and Ω
    Ω0  = deepcopy(Ωguess)
    Σ0  = deepcopy(Σguess)

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
        D = diagm(abs.(D).>1e-10).*diagm(D)
        U = (abs.(U).>1e-10).*U ;

        if maximum(abs.(imag.(D))) < 1e-6
            D = real.(D)
        end
        Λ       = U*max.(ω*eye - D,0.0)*U'; 
        Σ_p     = ω*SqRΣ*U/(max.(D,ω*eye))*U'*SqRΣ;
        
        Σ1      = A*Σ_p*A' + Q*Q'
        err         = norm(Σ1 - Σ0,2)/sqrt(n)
        
        Σ0      = w*Σ1 + (1-w)*Σ0
        
        SqRΣ    = sqrt(Σ0);
        invSqRΣ = pinv(SqRΣ);
        
        Ω0      = Ω_c + β*A'*invSqRΣ*(ω*eye - Λ)*invSqRΣ*A
        
        SqRΣ = (abs.(SqRΣ).>1e-10).*SqRΣ
        Ω0 = (abs.(Ω0).>1e-10).*Ω0

        iter += 1
    end
    
    inv_Σ1  = pinv(Σ1) ;
    Y = (eye - Σ_p*inv_Σ1)'*H ;
    Σ_z = H'*(Σ_p - Σ_p*inv_Σ1*Σ_p)*H ;
    K = Σ1*Y*pinv(Y'*Σ1*Y + Σ_z) ;

    

    return(Σ1,Σ_p,Λ,Ω0,Y,Σ_z,K)
end
