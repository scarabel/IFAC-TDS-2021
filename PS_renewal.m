function out = PS_renewal
% Copyright (c) 2021 Francesca Scarabel
% 
% MatCont system definition file of the PseudoSpectral Discretization of
% the nonlinear RE taken from Breda et al, EJQTDE 2016,
% x(t) = 0.5*gamma*int_a0^tau x(t-a)*exp(x(t-a))da
% for the integrated state B=int_0^t x(s)ds

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= []; %@userf;
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,loggamma,a0,tau,M) 

    % construction of Chebyshev nodes in [-tau,0] and differentiation matrix
    [~,Nodes,DD,BaryWeights] = cheb_delay(M,-tau,0);
    DM = DD(2:end,2:end);
    
    % compute derivative of the interpolating polynomial on the nodes
    der_state = DM*state;

    % construction of quadratures nodes and weights    
    [QuadWeights,QuadNodes,~,~] = cheb_delay(50,-tau,-a0);

    der = interpoly(QuadNodes,Nodes(2:end),der_state,Nodes(2:end)'.*BaryWeights(2:end));

    fM = 0.5*exp(loggamma) * (QuadWeights*(der.*exp(-der)));

    dydt = der_state - fM;
    
end
 
% --------------------------------------------------------------------------
function state_eq = init(M,xeq,tau)
% INPUT: M is the discretization parameter
%        xeq is the equilibrium value of the RE 
%        tau is the delay 
% OUTPUT state_eq is the initial vector for init_EP_EP

   [~,Nodes,~,~] = cheb_delay(M,-tau,0);
    state_eq = xeq*Nodes(2:end);
    
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

function pN = interpoly(Theta,Nodes,Values,Weights)
    % Computes the value of the interpolating polynomial (Nodes,Values) in theta,
    % using barycentric interpolation with Weights.

    N=length(Nodes)-1;
    L=length(Theta);
    pN=zeros(L,1);

    for ii=1:L
        diff=ones(1,N+1)*Theta(ii)-Nodes';
        if any(diff==0)
            pN(ii)=(diff==0)*Values;
        else        
        Coeff=Weights./diff;
        pN(ii) = (Coeff*Values)/sum(Coeff);
        end
    end

end
