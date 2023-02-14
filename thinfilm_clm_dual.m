% Description: This MATLAB program solves the thin-film equation 
%
%           h_t + (h^n h_xxx)_x = 0
%
% with nonzero contact angles and kinematic condition. The solution h(t,x) 
% is a function of space 'x' and time 't' and defined on a domain 
% omega=(x(1),x(n)) which depends on time. 
%
% This code extends to functionality of thinfilm_clm.m to the better 
% DUAL variational structure that handles dynamic contact angles and 
% is robust in partial/complete/dry wetting. This version allows general
% mobility exponents and different laws for dynamic contact angle.
%
% written by Dirk Peschka
% GNU GPL 2

clear all

% set model & computational paramters
SL     =    0.0;  % negative spreading coefficient at x=x-
SR     =    0.0;  % negative spreading coefficient at x=x+
g1     =    0.0;  % tangential gravity
g2     =    0.0;  % normal gravity

mobexp =    2.0;  % mobility exponent n in m(h)=h^n
L      =    2.0;  % initial domain size (0,L)

ppower =    2.0;  % power 2/alpha in dynamic contact angle law
mu0    =    1.0;  % prefactor in dynamic contact angle law
nt     =  10000;  % number of time steps
npoint =    512;  % number of vertices

% create element decomposition for FE method
x               = linspace(0,L,npoint)';% vertices
nelement        = npoint-1;             % no elements
nd(1:nelement,1)= 1:npoint-1;           % id left point of an element
nd(1:nelement,2)= 2:npoint;             % id right point of an element
local_mass_p1   = [1/3 1/6;1/6 1/3];    % mass matrix for reference [0,1]

% * create & remember initial data
h  = (1-(x-1).^2);
t  = 0;
dt = 1d-3;
integration_weights_dual

for it=1:nt 
    % plot numerical solution every 1000-th iteration
    if mod(it,100)==1
        fprintf('it=%i time=%e vol=%f \n',it,t,sum(diff(x).*(h(1:end-1)+h(2:end))/2))
        plot(x,h,'b-','LineWidth',2);
        xlim([-8 8])
        ylim([0 2])
        drawnow
    end
    
    % construct system matrices
    build_FE_matrices_dual % script: matrices A,S,M,Dx for FEM
    build_ALE_matrix_dual         % script: matrices for ALE decomposition
    
    % Mbd = sparse([1 2],[1 ndof],[1/abs(dh(1)) 1/abs(dh(end))]);
    % is the 1D boundary mass matrix
    Mbd = sparse([1 2],[1 ndof],[1 1]);
    ZZ  = sparse(2,ndof);
    ZZ1 = sparse(2,2);
    Sbd = mu0*[abs(dh(1))^ppower 0;0 abs(dh(end))^ppower];
    
    % FE problem: build right-hand-side rhs 
    rhs=[S*h+M*(2*g2*h-g1*x);zeros(ndof+2,1);];
    rhs(1)=rhs(1)+(SL+(dh(  1)^2)/2)/(dh(  1));
    rhs(ndof)=rhs(ndof)-(SR+(dh(end)^2)/2)/(dh(end));
    
    % construction of gradient flow saddle point structure with 
    % bulk mobility in Sw and contact line mobility in Sbd
    % the term "-dt*S" emerges from semi-implicit treatment of 
    % the fourth-order term
    A = [-dt*S  M Mbd';...
        M Sw  ZZ';...
        Mbd ZZ Sbd];
    
    u = A\rhs; % solve for u=(hdot,pi)^t
    
    % perform ALE decomposition & update solution
    U = (I-P)\u(1:ndof); % select only u, forget p
    h = h + dt*I*U;      % update h
    x = x + dt*X*U;      % update x
    t = t + dt;  
    
    % problem specific postprocessing 
    % and timestepping needs to be added here
end