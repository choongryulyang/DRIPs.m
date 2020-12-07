%% Afrouzi and Yang (2019)
% This code returns the irfs of a Drip either with the steady state
% information structure or in the transition path if specified.
%
%% SYNTAX
%    irfs(p, kwargs...) -> Path
%
%   where p is the output of the Drip.m function
%
%% OPTIONAL ARGUMENTS
%   default values are set unless specified otherwise by the user 
%   trip   = none  : irfs under the specified transition path as the output of Trip.m function
%                    if no transition path is specified, returns irfs under the steady state 
%                    information structure by default.
%   length = 40    : length of irfs
%
%% OUTPUTS
%   returns a structure `out` that contains the following fields:
%   out.x    : this is the impulse responses of x itself to structural
%              shocks and has dimensions n*k*T where `n` is the dimension
%              of x, `k` is the numnber of structural shocks and `T` is the
%              length of irfs.
%   out.x_hat: this is the impulse responses of beliefs to shocks and has 
%              dimension n*k*T where `n` is the dimension of the state
%              vector x, `k` is the number of structural shocks, and `T` is
%              the length of irfs. 
%   out.a    : this is the impulse respnses of actions to shocks and has
%              dimension m*k*T where `m` is the number of actions, `k` is
%              the number of structural shocks and `T` is the length of
%              irfs. For instance a(i,j,:) is the irf of i'th action to
%              j'th structural shock.
% 
%% EXAMPLES:
%  % irfs with the steady state information structure
%       >> p     = Drip(omega,beta,A,Q,H) % solve a Drip
%       >> pirfs = irfs(p) 
%   
%  % irfs with the transition dynamics in information structure
%       >> p_trip = Trip(p,Sigma0) % solve for a transition path given Sigma0
%       >> pirfs = irfs(p,'trip',p_trip)
%
%  % adjusting the length of irfs 
%       >> pirfs = irfs(p,'length',20) % returns irfs for 20 periods 
%                                      % (default is 40)

function out = irfs(p,varargin)
    args = inputParser; 
    % parse optional inputs
    addOptional(args,'trip','none');
    addOptional(args,'length',40);

    parse(args,varargin{:});
    
    [n,m]      = size(p.H) ;
    [~,k]      = size(p.Q) ;   
    Ik         = eye(k);

    if any(strcmp(args.UsingDefaults,'trip'))
        out = dripirfs(p,args.Results.length);
    else 
        out = tripirfs(p,args.Results.trip,args.Results.length);
    end 

    function out = dripirfs(p,length) 
        % initialize 
        out.T     = length;
        out.x     = zeros(n,k,out.T);
        out.x_hat = out.x;
        out.a     = zeros(m,k,out.T);

        for kk = 1:1:k
            e_k = Ik(:,kk);
            for ii = 1:1:out.T
                if ii == 1
                    out.x(:,kk,ii) = p.Q*e_k;
                    out.x_hat(:,kk,ii) = (p.ss.K*p.ss.Y')*(out.x(:,kk,ii));
                else
                    out.x(:,kk,ii) = p.A*out.x(:,kk,ii-1);
                    out.x_hat(:,kk,ii) = p.A*out.x_hat(:,kk,ii-1) + ...
                        (p.ss.K*p.ss.Y')*(out.x(:,kk,ii)-p.A*out.x_hat(:,kk,ii-1));
                end
                out.a(:,kk,ii) = p.H'*out.x_hat(:,kk,ii);
            end
        end 
    end

    function out = tripirfs(p,pt,length) 
        %initialize
        out.T     = length;
        L         = size(pt.Ds,2);
        
        out.x     = zeros(n,k,out.T);
        out.x_hat = out.x;
        out.a     = zeros(m,k,out.T);
        In        = eye(n);

        for kk = 1:1:k
            e_k = Ik(:,kk);
            for ii = 1:out.T
                if ii==1
                    out.x(:,kk,ii)     = p.Q*e_k;
                    out.x_hat(:,kk,ii) = (In-pt.Sigma_ps(:,:,ii)* ...
                        pinv(pt.Sigma_1s(:,:,ii)))*(out.x(:,kk,ii));
                elseif ii <= L
                    out.x(:,kk,ii)     =p.A*out.x(:,kk,ii-1);
                    out.x_hat(:,kk,ii) =p.A*out.x_hat(:,kk,ii-1)+ ...
                        (In-pt.Sigma_ps(:,:,ii)*pinv(pt.Sigma_1s(:,:,ii)))* ...
                        (out.x(:,kk,ii)-p.A*out.x_hat(:,kk,ii-1));
                else
                    out.x(:,kk,ii)     =p.A*out.x(:,kk,ii-1);
                    out.x_hat(:,kk,ii) =p.A*out.x_hat(:,kk,ii-1)+ ...
                        (In-pt.Sigma_ps(:,:,end)*pinv(pt.Sigma_1s(:,:,end)))* ...
                        (out.x(:,kk,ii)-p.A*out.x_hat(:,kk,ii-1));
                end
                out.a(:,kk,ii)  = p.H'*out.x_hat(:,kk,ii);
            end
        end
    end  
end

