% Description: This MATLAB program solves the thin-film equation
%           h_t + (h^2 h_xxx)_x = 0
% with nonzero contact angles and kinematic condition. The solution h(t,x)
% is a function of space 'x' and time 't' and defined on a domain
% omega=(x(1),x(n)) which depends on time. Details can be found in the
% corresponding paper xxx in Journal of Computational Physics (submitted)
%
% Added contact line model
%
% written by Dirk Peschka
% GNU GPL 2

clear all

% set model & computational paramters
L     =  1.0;  % initial domain size (0,L)
T     =  1*4.0;  % final time
SL    =  1.0;  % negative spreading coefficient at x=x-
SR    =  1.0;  % negative spreading coefficient at x=x+
g1    = 20.0;  % tangential gravity
g2    =  0.0;  % normal gravity
nt    =  10000;  % number of time steps
npoint=  1*500;  % number of vertices
beta  =  3/4;
% create element decomposition for FE method
x               =linspace(0,L,npoint)';% vertices
nelement        =npoint-1;             % no elements
nd(1:nelement,1)=1:npoint-1;           % id left point of an element
nd(1:nelement,2)=2:npoint;             % id right point of an element
local_mass_p1   =[1/3 1/6;1/6 1/3];    % mass matrix for reference [0,1]

% * create & remember initial data
h  = x.*(L-x);

hi = h;xi = x;
t  = 0;dt = T/nt;

% initial data
%hold off
%plot(xi,hi,'r--','LineWidth',2);
%hold on
hold off
m=0;
for it=1:nt

    % construct system matrices
    build_FE_matrices % script: matrices A,S,M,Dx for FEM
    build_ALE_matrix  % script: matrices for ALE decomposition
    % plot numerical solution
    if mod(it,10)==1
        plot(xi,hi,'r--',x,h,'b-','LineWidth',0.5);
        xlim([0 5])
        ylim([0 0.5])
        view(2)
        box on
        drawnow
    end

    % contact line dissipation
    DD = beta*[1 0;0 1];
    
    % different constraints
    C1 = [1 zeros(1,ndof-1)];
    C2 = [zeros(1,ndof-1) 1];
    C3 = [dh(1)   0];
    C4 = [0 dh(end)];
    C5 = ones(1,ndof)*M;
    
    % zeros of different sizes
    Z1 = zeros(1,ndof);
    Z2 = zeros(2,ndof);
    Z3 = [0 0];
    
    % u = (hdot,pi,xdot,lambda)
    D  = [dt*S  Z Z2'  M C1' C2' Z1';
             Z Sw Z2' Sw Z1' Z1' C5';
            Z2 Z2 DD  Z2 C3' C4' Z3';
             M Sw Z2'  Z Z1' Z1' Z1';
            C1 Z1 C3  Z1  0   0   0 ;
            C2 Z1 C4  Z1  0   0   0 ;
            Z1 C5 Z3  Z1  0   0   0];

    % FE problem: build right-hand-side rhs & solve
    rhs=[-S*h-M*(2*g2*h-g1*x);zeros(ndof,1)];
    rhs=[rhs;+(SL+(dh(  1)^2)/2)];
    rhs=[rhs;-(SR+(dh(end)^2)/2)];
    rhs=[rhs;zeros(ndof+2,1)];
    rhs=[rhs;0];
    
    % solve
    u = D\rhs;
    
    % perform ALE decomposition & update solution
    U = (I-P)\u(1:ndof); % select only u, forget p
    h = h + dt*I*U;      % update h
    x = x + dt*X*U;      % update x
    t = t + dt;


end