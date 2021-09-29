% copyright(c) 2021 Francesca Scarabel
%
% MC_logistic
% MatCont command line instructions for the numerical bifurcation analysis
% of the delayed logistic equation defined in PS_logistic.m

clear;
close all;

% Initial parameter values
r = 0.1; 
M = 10;
par = [r,M]';

ap1 = 1; % index of the continuation parameter in the vector par
handles = feval(@PS_logistic);

%% Equilibrium continuation from initial point

% set options for MatCont continuation
opt=contset;
opt=contset(opt,'MaxNumPoints',200);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);

init_eq = ones(M+1,1); % initial estimate of equilibrium point
[x0,v0] = init_EP_EP(@PS_logistic,init_eq,par,ap1);
[xe,ve,se,he,fe] = cont(@equilibrium,x0,v0,opt); 
% last row of xe contains the parameter values
% rows 1:(M+1) contain the approximation of the equilibrium values
xe(end,end) % last computed parameter value

% global cds;
% [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);

% Plot equilibrium curve
figure(1); clf;
cpl(xe,ve,se,[(M+1)+1 1]);
title('Bifurcation diagram w.r.t. $r$','Interpreter','latex');
xlabel('$r$','interpreter','latex');

%% Print eigenvalues corresponding to a given point
figure
j = se(2).index;
fe(:,j)
subplot(1,2,1)
plot(real(fe(:,j)),imag(fe(:,j)),'o')
xline(0,'k--'); yline(0,'k--');
axis([-30 5 -25 25])
title(['Eigenvalues at $r=$',num2str(xe(end,j))],'interpreter','latex')
j = 50;
fe(:,j)
subplot(1,2,2)
plot(real(fe(:,j)),imag(fe(:,j)),'o')
xline(0,'k--'); yline(0,'k--');
axis([-30 5 -25 25])
title(['Eigenvalues at $r=$',num2str(xe(end,j))],'interpreter','latex')

%% Save data for detected Hopf point and limit cycle continuation
H_index = se(2).index; % index of Hopf point in equilibrium branch
H = xe(1:M+1,H_index); % state vector at Hopf point
par(ap1) = xe(end,H_index);

% Limit cycle continuation from H
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'MaxStepsize',1);
ntst = 20; % number of intervals in one period
ncol = 2; % degree of collocation polynomial

[x0,v0] = init_H_LC(@PS_logistic,H,par,ap1,1e-6,ntst,ncol);
[xlc,vlc,slc,hlc,flc] = cont(@limitcycle,x0,v0,opt); 
xlc(end,end) % last computed parameter value
% [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)

% Plot max and min of periodic solutions
upperbound = max(xlc(1:(M+1):((ntst*ncol+1)*(M+1)),:));
lowerbound = min(xlc(1:(M+1):((ntst*ncol+1)*(M+1)),:));
figure(1); hold on
plot(xlc(end,:),upperbound,'b',xlc(end,:),lowerbound,'b');
axis([0 3 0 5])

% Plot profile of orbit corresponding to index j
% last row in xlc contains the parameter value
% row 'end-1' in xlc contains the approximated period
% rows 1:ntst+1 in flc contain the points of the solution approximated in one period (normalised to 1)
j=30;
figure
mesh= xlc(end-1,j)*flc(1:ntst+1,j);
profile= xlc(1:(M+1)*ncol:((ntst*ncol+1)*(M+1)),j);
plot(mesh,profile);
xlim([0,xlc(end-1,j)]); 
title(['Periodic solution corresponding to $r=$',num2str(xlc(end,j),3)],'Interpreter','latex');
