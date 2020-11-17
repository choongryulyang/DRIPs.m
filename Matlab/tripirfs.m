%% Afrouzi and Yang (2019)
% This code returns the irfs of the a Drip under the transition dynamics of 
% the information structure from an initial belief x_0 ~ N(0,Sigma0)
% see the Julia package for transition from non-zero mean for prior beliefs
%% SYNTAX
%     tripirfs(p, kwargs...) -> Path
%
%   where p is the output of the Trip.m function
%
%% OPTIONAL ARGUMENTS
%   default values are set unless specified otherwise by the user 
%   length = 40    : length of irfs 

function irfs = tripirfs(pt,varargin) 

    % parse optional inputs
    args = inputParser;
    addOptional(args,'length',40);

    parse(args,varargin{:});

    % initialize 
    irfs.T     = args.Results.length;
    p          = pt.p;
    [n,m]      = size(p.H) ;
    [~,k]      = size(p.Q) ;
    L          = size(pt.Ds,2);
    
    irfs.x     = zeros(n,k,irfs.T);
    irfs.x_hat = irfs.x;
    irfs.a     = zeros(m,k,irfs.T);
    In         = eye(n);
    Ik         = eye(k);

    for kk = 1:1:k
        e_k = Ik(:,k);
        for ii = 1:irfs.T
            if ii==1
                irfs.x(:,kk,ii)     = p.Q*e_k;
                irfs.x_hat(:,kk,ii) = (In-pt.Sigma_ps(:,:,ii)* ...
                    pinv(pt.Sigma_1s(:,:,ii)))*(irfs.x(:,kk,ii));
            elseif ii <= L
                irfs.x(:,kk,ii)     =p.A*irfs.x(:,kk,ii-1);
                irfs.x_hat(:,kk,ii) =p.A*irfs.x_hat(:,kk,ii-1)+ ...
                    (In-pt.Sigma_ps(:,:,ii)*pinv(pt.Sigma_1s(:,:,ii)))* ...
                    (irfs.x(:,kk,ii)-p.A*irfs.x_hat(:,kk,ii-1));
            else
                irfs.x(:,kk,ii)     =p.A*irfs.x(:,kk,ii-1);
                irfs.x_hat(:,kk,ii) =p.A*irfs.x_hat(:,kk,ii-1)+ ...
                    (In-pt.Sigma_ps(:,:,end)*pinv(pt.Sigma_1s(:,:,end)))* ...
                    (irfs.x(:,kk,ii)-p.A*irfs.x_hat(:,kk,ii-1));
            end
            irfs.a(:,kk,ii)  = p.H'*irfs.x_hat(:,kk,ii);
        end
    end 
end 
