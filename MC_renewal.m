% Copyright (c) 2021 Francesca Scarabel
% 
% MC_renewal.m
% MatCont command line instructions for the numerical bifurcation analysis
% of the renewal equation for cannibalism defined in PS_renewal.m

clear;
close all;

% Initial parameter values
tau = 3;
abar = 1;
loggamma = -1;
M = 10;
par = [loggamma,abar,tau,M]';
ap1 = 1; % index of the continuation parameter in the vector par
handles = feval(@PS_renewal);

%% Equilibrium continuation from initial point [xeq;yeq]

% set options
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',50);

% Approximated equilibrium corresponding to par
xeq = 0;
state_eq = feval(handles{1},M,xeq,tau); % initializes equilibrium vector

[x0,v0] = init_EP_EP(@PS_renewal,state_eq,par,ap1);
[xe,ve,se,he,fe] = cont(@equilibrium,x0,v0,opt); 
xe(end,end)
% global cds;
% [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);

figure(10); clf;
cpl(xe,ve,se,[M+1 M]);
xlabel('$\log(\gamma)$', 'Interpreter','latex');
hold on;

%% Equilibrium continuation from BP

% Detection of branching point
sBP = se(2);
BP_index = se(2).index;
BP = xe(1:M,BP_index);
par(ap1) = xe(end,BP_index);

% set options
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',500);

[x0,v0] = init_BP_EP(@PS_renewal,BP,par,sBP,0.1);
[xe,ve,se,he,fe] = cont(@equilibrium,x0,v0,opt); xe(end,end)
% [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)

% % Plot
% figure(1)
% cpl(xe,ve,se,[M+1 M]);
%
% Plot of bifurcation diagram of first component
[~,Nodes,DD,BaryWeights] = cheb_delay(M,-tau,0);
DM = DD(2:end,2:end);
[QuadWeights,QuadNodes,~,~]=cheb_delay(50,-tau,-abar);

for index_sol = 1:size(xe,2)
    BB = DM*xe(1:end-1,index_sol);
    der = interpoly(QuadNodes,Nodes(2:end),BB,Nodes(2:end)'.*BaryWeights(2:end));
    FM = 0.5*exp(xe(end,index_sol))*QuadWeights*(der.*exp(-der));
    b0(index_sol) = FM;
end

figure(10)
plot(xe(end,:),b0,'b');
hold on

for ii=2:length(se)-1
    index=se(ii).index;
    plot(xe(end,index),b0(index),'.r');
end

%% H continuation in two parameters

% Detection of Hopf bifurcation point
H_index = se(2).index;
H = xe(1:M,H_index);
par(ap1) = xe(end,H_index);

ap2 = 3; % index of second parameter for two-parameter continuation
TOL = 1e-3;
% set options
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'Singularities',0);
opt=contset(opt,'Eigenvalues',0);

% forward
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',1);

[x0,v0]=init_H_H(@PS_renewal,H,par,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(M+1,end)
% [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)

% Plot
figure
cpl(xh,vh,sh,[M+1 M+2]); hold on
axis([0 10 0 3])
xlabel('$\log(\gamma)$','Interpreter','latex'); ylabel('$\tau$','Interpreter','latex');
title('Stability regions (Hopf bifurcation curve)')

% % backward
% opt=contset(opt,'Backward',1);
% [x0,v0]=init_H_H(@PS_renewal,H,parH,[ap1 ap2]);
% [xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
% [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
% cpl(xhb,vhb,shb,[MM+1 MM+2]); hold on

%% Limit cycle continuation from H
ntst=20; % number of interval
ncol=2; % degree of polynomial

TOL=1e-3;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

% set options
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',1);
opt=contset(opt,'InitStepsize',0.1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',0);
opt=contset(opt,'MaxNumPoints',50);

[x0,v0]=init_H_LC(@PS_renewal,H,par,ap1,1e-3,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
% [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)

%% Plot max and min periodic solutions

Per_Solutions = zeros(ntst*ncol+1,size(xlc,2));
for ind_persol=1:size(xlc,2)
    for ind_mesh=1:ntst*ncol+1
        BB = DM*xlc((ind_mesh-1)*M+1:(ind_mesh-1)*M+M,ind_persol);
        der = interpoly(QuadNodes,Nodes(2:end),BB,Nodes(2:end)'.*BaryWeights(2:end));
        b0_per = 0.5*exp(xlc(end,ind_persol))*QuadWeights*(der.*exp(-der));
        Per_Solutions(ind_mesh,ind_persol) = b0_per;
    end
end

upperbound=max(Per_Solutions);
lowerbound=min(Per_Solutions);

figure(10); hold on
plot(xlc(end,:),upperbound,'b',xlc(end,:),lowerbound,'b');

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'or',xlc(end,index),lowerbound(index),'or');
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