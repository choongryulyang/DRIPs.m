%% Afrouzi and Yang (2019)
% This code returns the irfs of the a Drip under the steady state 
% information structure
%
%% SYNTAX
%     dripirfs(p, kwargs...) -> Path
%
%   where p is the output of the Drip.m function
%
%% OPTIONAL ARGUMENTS
%   default values are set unless specified otherwise by the user 
%   length = 40    : length of irfs 

function irfs = dripirfs(p,varargin) 

    % parse optional inputs
    args = inputParser;
    addOptional(args,'length',40);

    parse(args,varargin{:});

    % initialize 
    irfs.T     = args.Results.length;
    [n,m]      = size(p.H) ;
    [~,k]      = size(p.Q) ;
    irfs.x     = zeros(n,k,irfs.T);
    irfs.x_hat = irfs.x;
    irfs.a     = zeros(m,k,irfs.T);
    Ik         = eye(k);

    for kk = 1:1:k
        e_k = Ik(:,k);
        for ii = 1:1:irfs.T
            if ii == 1
                irfs.x(:,kk,ii) = p.Q*e_k;
                irfs.x_hat(:,kk,ii) = (p.ss.K*p.ss.Y')*(irfs.x(:,kk,ii));
            else
                irfs.x(:,kk,ii) = p.A*irfs.x(:,kk,ii-1);
                irfs.x_hat(:,kk,ii) = p.A*irfs.x_hat(:,kk,ii-1) + ...
                    (p.ss.K*p.ss.Y')*(irfs.x(:,kk,ii)-p.A*irfs.x_hat(:,kk,ii-1));
            end
            irfs.a(:,kk,ii) = p.H'*irfs.x_hat(:,kk,ii);
        end
    end 
end 
