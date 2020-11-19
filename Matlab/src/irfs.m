%% Afrouzi and Yang (2019)
% This code returns the irfs of the a Drip 
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
        L          = size(pt.Ds,2);
        
        out.x     = zeros(n,k,out.T);
        out.x_hat = out.x;
        out.a     = zeros(m,k,out.T);
        In         = eye(n);

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

