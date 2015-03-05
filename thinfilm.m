% Description: This MATLAB program solves the thin-film equation 
%           h_t + (h^2 h_xxx)_x = 0
% with nonzero contact angles and kinematic condition. The solution h(t,x) 
% is a function of space 'x' and time 't' and defined on a domain 
% omega=(x(1),x(n)) which depends on time. Details can be found in the 
% corresponding paper xxx in Journal of Computational Physics (submitted)
%
% written by Dirk Peschka
% GNU GPL 2

clear all

% set model & computational paramters
L     = 1.0;  % initial domain size (0,L)
T     = 0.2;  % final time
SL    = 1.0;  % negative spreading coefficient at x=x-
SR    = 1.0;  % negative spreading coefficient at x=x+
g1    = 0.0;  % tangential gravity
g2    = 0.0;  % normal gravity
nt    = 100;  % number of time steps
npoint= 100;  % number of vertices

% create element decomposition for FE method
x               =linspace(0,L,npoint)';% vertices
nelement        =npoint-1;             % no elements
nd(1:nelement,1)=1:npoint-1;           % id left point of an element
nd(1:nelement,2)=2:npoint;             % id right point of an element
local_mass_p1   =[1/3 1/6;1/6 1/3];    % mass matrix for reference [0,1]

% * create & remember initial data
h  = L/2-abs(L/2-x); 

hi = h;xi = x;
t  = 0;dt = T/nt;

% initial data
hold off
plot(xi,hi,'r--','LineWidth',2);
xlabel('x','FontSize',22);ylabel('h','FontSize',22);
hold on
for it=1:nt       
    % construct system matrices
    build_FE_matrices % script: matrices A,S,M,Dx for FEM
    build_ALE_matrix  % script: matrices for ALE decomposition
    
    % FE problem: build right-hand-side rhs & solve
    rhs=[zeros(npoint,1);S*h+M*(2*g2*h-g1*x)];
    rhs(ndof+1)=rhs(ndof+1)+(SL+(dh(  1)^2)/2)/dh(  1);
    rhs(2*ndof)=rhs(2*ndof)-(SR+(dh(end)^2)/2)/dh(end);
    u = A\rhs; % solve for u=(hdot,pi)^t
    
    % perform ALE decomposition & update solution
    U = (I-P)\u(1:ndof); % select only u, forget p
    h = h + dt*I*U;      % update h
    x = x + dt*X*U;      % update x
    
    % plot numerical solution
    if mod(it,10)==1
        plot(x,h,'b-','LineWidth',2);
        drawnow
    end
end
legend('initial','t > 0')
set(gca,'FontSize',22)
hold off