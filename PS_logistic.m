function out = PS_logistic
% Matcont system definition file of the PseudoSpectral Discretization of
% the delayed logistic equation
% y'(t) = r*y(t)*[1-y(t-tau)]

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= []; 
out{11}= []; %@userf1;
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,r,M) 
% OUTPUT: evaluation of the rhs of the approximating ODE system
% INPUT: 
% state - vector of variables, dimension (M+1)
% M - discretization index

    % construction of Chebyshev nodes in [-1,0] and differentiation matrix
    [QuadWeights,Nodes,DD,BaryWeights] = cheb_delay(M,-1,0);

    fM = r*state(1).*(1-state(end));

    dydt = [ fM; DD(2:end,:)*state ];

end
 
 
% --------------------------------------------------------------------------
function state_eq = init(M,yeq)
% INPUT: M is the discretization parameter
%        yeq is the equilibrium value of the DDE 
% OUTPUT state_eq is the initial vector for init_EP_EP
    
    state_eq = yeq*ones(M+1,1);
    
end

% ------------
%%% AUXILIARY FUNCTIONS

function [w,x,D,q]=cheb_delay(N,a,b)
    % see Trefethen, Spectral Methods in Matlab, SIAM, 2000
    % Output:
    % x - N+1 Chebyshev nodes on [a,b] (x_0=b, x_N=a),
    % w - weights of the quadrature formula in [a,b],
    % D - differentiation matrix
    % q - row vector of the barycentric weights

    if N==0
        x=1;
        D=0;
        return
    end
    p=pi*(0:N)'/N;
    x=((b-a)*cos(p)+b+a)/2;

    c=[2;ones(N-1,1);2].*(-1).^(0:N)';
    X=x(:,ones(1,N+1)); %X=repmat(x,1,N+1);
    dX=X-X';
    D=(c*(1./c)')./(dX+(eye(N+1)));
    D=D-diag(sum(D,2)); %D=D-diag(sum(D'));

    % Quadrature weights
    w=zeros(1,N+1);
    ii=2:N;
    v=ones(N-1,1);
    if mod(N,2)==0
        w(1)=1/(N^2-1);
        w(N+1)=w(1);
        for k=1:N/2-1
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
        v=v-cos(N*p(ii))/(N^2-1);
    else
        w(1)=1/N^2;
        w(N+1)=w(1);
        for k=1:(N-1)/2
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
    end
    w(ii)=2*v/N;
    w=w*abs(b-a)/2;

    % Barycentric weights
    q=1./prod(dX'+eye(N+1)); %q=1./prod(dX'+eye(N+1)); % row vector of the barycentric weights
end