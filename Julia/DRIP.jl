## A Software Package for Solving Dynamic Rational Inattention Problems (D.R.I.P.)
## Copyright © 2020 Hassan Afrouzi, Choongryul Yang and Miguel Acosta
## When using this package for your work, cite Afrouzi and Yang (2019) for reference.
## Julia code written by Miguel Acosta / Modified by Hassan Afrouzi

using LinearAlgebra

## A General Structure for D.R.I.P.
struct drip 
    ω; β; A; Q; H;                    # primitives
    K; Y; Σ_z; Σ_p; Σ_1; Ω;           # solution
end

## A General Structure for Transition dynamics of Rational Inattention Problems (T.R.I.P.)
struct trip 
    P::drip;                # problem and solution in steady state
    T::Int;                 # length of T.R.I.P.
    Σ_1s; Σ_ps; Ωs;         # priors, posteriors and benefit matrices 
    Ds;                     # eigenvalues and eigenvectors of Σ_t^(0.5)Ω_tΣ_t^(0.5)
end 

## A Structure for the impulse responses of D.R.I.P. for the steady state posterior
struct dripirfs
    T     :: Int            # length of IRFs
    x     :: Array{Float64} # IRFs of the fundamental shocks
    x_hat :: Array{Float64} # IRFs of beliefs
    a     :: Array{Float64} # IRFs of actions
end

## solve_drip: this function solves for the steady state of of the D.R.I.P.
## Inputs
# ω      : Cost of information
# β      : Discount factor
# A      : Transition matrix: x=Ax+Qu
# Q      : Std. Dev. matrix: x=Ax+Qu
# H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)

## Outputs: a drip structure with 
# Σ_1    : Steady-state Prior uncertainty
# Σ_p    : Steady-state posterior uncertainty
# Λ      : Shadow matrix on the no-forgetting constraint -- Λ*(Σ-Σ_p) = 0
# Ω      : Dynamic benefit matrix
# Y      : Weight vector for evolution of actions
# Σ_z    : Covariance matrix of the rational inattention error
# K      : Kalman gain matrix

function solve_drip(ω,β,A,Q,H;              # primitives of the D.R.I.P.
                    fcap::Bool = false,     # optional: if true then solves the problem with fixed capacity κ = ω.
                    Ω0         = H*H',      # optional: initial guess for steady state information matrix
                    Σ0         = A*A'+Q*Q', # optional: initial guess for steady state prior
                    w          = 1,         # optional: updating weight in iteration
                    tol        = 1e-4,      # optional: tolerance level for convergence
                    maxit      = 10000)     # optional: maximum number of iterations

    ## initialize
    (n,m) = length(size(H)) == 2 ? size(H) : (size(H,1),1);  
    # n: dimension of state, m: number of actions
    eye   = Matrix{Float64}(I,n,n);
    err   = 1
    iter  = 0
    SqRΣ  = sqrt(Σ0)
    Ω_c   = H*H';
    Σ1    = Matrix{Float64}(I,n,n)
    Σ_p   = Matrix{Float64}(I,n,n)
    Λ     = Matrix{Float64}(I,n,n)
    κ     = ω
    # iterate 
    while (err > tol) & (iter < maxit)
        D, U    = eigen(SqRΣ*Ω0*SqRΣ);
        D       = diagm(getreal(D));
        U       = getreal(U);
        
        if fcap == true 
            ω = (2^(2*κ)/det(max.(ω*eye,D)))^(-1/n);
        end
        Λ       = U*(max.(ω*eye-D,0.0))*U';
        Σ_p     = getreal(ω*SqRΣ*U/(max.(D,ω*eye))*U'*SqRΣ);

        Σ1      = A*Σ_p*A' + Q*Q'
        err     = norm(Σ1 - Σ0,2)/norm(Σ0,2);

        Σ0      = w*Σ1 + (1-w)*Σ0

        SqRΣ    = sqrt(Σ0);
        invSqRΣ = pinv(SqRΣ);

        Ω0      = w*(Ω_c + β*A'*invSqRΣ*(ω*eye - Λ)*invSqRΣ*A)+(1-w)*Ω0;
        Ω0      = (abs.(Ω0).>1e-10).*Ω0;

        iter   += 1
    end

    if iter == maxit
        print("RI Code hit maxit\n")
    end

    Σ1     = (abs.(Σ1).>1e-10).*Σ1 + diagm(abs.(diag(Σ1)).<= 1e-10)*1e-8
    inv_Σ1 = pinv(Σ1) ;
    Y      = (eye - Σ_p*inv_Σ1)'*H ;
    Σ_z    = H'*(Σ_p - Σ_p*inv_Σ1*Σ_p)*H ;
    K      = Σ1*Y*pinv(Y'*Σ1*Y + Σ_z) ;

    return(drip(ω,β,A,Q,H,K,Y,Σ_z,Σ_p,Σ1,Ω0))
end

## solve_trip: this function solves for the transition dynamics given an initial prior

function solve_trip(P::drip,         # D.R.I.P.
                    Σ0;              # initial prior matrix
                    T     = 100,     # optional: time until convergence to steady state
                    tol   = 1e-4,    # optional: tolerance for convergence
                    maxit = 1000     # optional: max iterations
                    )
    ## Initialize 
    Ωs  = repeat(P.Ω, inner = [1,1,T]);
    Ωsp = Ωs;

    Σ_1s  = repeat(P.Σ_1, inner = [1,1,T]); Σ_1s[:,:,1] = Σ0;
    Σ_1sp = Σ_1s;
    (n,m) = length(size(P.H)) == 2 ? size(P.H) : (size(H,1),1);  
    # n: dimension of state, m: number of actions

    Σ_ps= repeat(P.Σ_p, inner = [1,1,T]); # initialize posteriors
    Ds  = zeros(n,T);                          # initialize eigenvalues 

    iter = 0;
    err  = 1;
    eye  = Matrix(I,n,n);
    while (err > tol) & (iter <= maxit)
        ## Given Ωs, find Sigmas using the law of motion for priors 
        for i in 1:1:T-1
            SqSigma        = real.(sqrt(Σ_1s[:,:,i]));
            D, U           = eigen(SqSigma*Ωs[:,:,i]*SqSigma);
            Ds[:,i]        = real.(D);
            U              = real.(U);
            D              = diagm(Ds[:,i]);
            Σ_ps[:,:,i]    = P.ω*SqSigma*U/(max.(D,P.ω*eye))*U'*SqSigma;
            Σ_1sp[:,:,i+1] = P.Q*P.Q'+P.A*Σ_ps[:,:,i]*P.A';
        end
        # Given Σ_1s, Find Omegas using the Euler equation
        for i = T-1:-1:1
            SqSigma      = sqrt(Σ_1sp[:,:,i+1]);
            invSqSigma   = real.(pinv(SqSigma));
            D, U         = eigen(SqSigma*Ωs[:,:,i+1]*SqSigma);
            D            = diagm(getreal(D));
            U            = real.(U);
            Ωsp[:,:,i]   = P.H*P.H'+P.β*P.A'*invSqSigma*U*min.(D,P.ω*eye)*U'*invSqSigma*P.A;
        end
        err    = norm(Σ_1sp-Σ_1s)/norm(Σ_1s);
        Σ_1s   = real.(Σ_1sp);
        Ωs     = real.(Ωsp);
        iter  += 1;
    end
    # Throw error if T was too short for convergence to steady state
    con_err = norm(P.Q*P.Q'+P.A*Σ_ps[:,:,end-1]*P.A'-P.Σ_1)/norm(P.Σ_1);
    if (con_err  > tol)
        error("T was too short for convergence. Try larger T.");
    end 
    # Finally, store the eigenvalues of steady state
    SqSigma      = real.(sqrt(Σ_1s[:,:,end]));
    D, U         = eigen(SqSigma*Ωs[:,:,end]*SqSigma);
    Ds[:,end]    = real.(D);
    # return the trip structure 
    return(trip(P,T,Σ_1s,Σ_ps,Ωs,Ds))
end

########## Aux. Functions ############

function getreal(M)
    #if maximum(abs.(imag.(M))) > 1e-6
    #    print("Warning: Matrix decomposition returned complex numbers larger than 1e-6")
    #end
    return(real.(M))
end

function infinitesum(func; tol = 1e-6,maxit = 1000,start=0)
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

## calculate the amount of information acquired in bits/nats for a drip 
function capacity(P::drip;      # drip structure
                  unit = "bit"  # optional: unit of capacity (bit or nat). 
                  )
    if unit == "bit"
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p))/log(2); #returns capacity in bits
    elseif unit == "nat" 
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p));        #returns capacity in nats 
    else 
        println("Invalid input for unit! Capacity is reported in bits.")
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p))/log(2);
    end
    return(κ)
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

function dripirfs(Pt::trip,T::Int)
    Ss    = Pt.P; # steady state drip
    (n,m) = length(size(Ss.H)) == 2 ? size(Ss.H) : (size(Ss.H,1),1) # dimensions of state and actions
    (_,k) = length(size(Ss.Q)) == 2 ? size(Ss.Q) : (size(Ss.Q,1),1) # number of structural shocks
    eye   = Matrix(I,n,n);

    x     = zeros(n,k,T); # initialize IRFs of state 
    x_hat = zeros(n,k,T); # initialize IRFs of beliefs 
    a     = zeros(m,k,T); # initialize IRFs of actions
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1 
                x[:,kk,ii]     = Ss.Q*e_k;
                x_hat[:,kk,ii] = (eye-Pt.Σ_ps[:,:,ii]*pinv(Pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]);
            else 
                x[:,kk,ii]     =Ss.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(eye-Pt.Σ_ps[:,:,ii]*pinv(Pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= Ss.H'*x_hat[:,kk,ii];
        end
    end
    return(dripirfs(T,x,x_hat,a))
end