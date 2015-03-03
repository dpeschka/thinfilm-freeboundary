% Solves the thin-film equation
%           h_t + (h^2 h_xxx)_x = 0
% with nonzero contact angles and kinematic condition.
clear all

% * set model & computational paramters
L     =  1.0;  % initial domain size (0,L)
T     =  1.0;  % final time
SL    =  1.0;  % negative spreading coefficient at x=x-
SR    =  1.0;  % negative spreading coefficient at x=x+
c1    =  0.0;  % normal gravity
c2    =  0.0;  % tangential gravity
nt    =  200;  % number of time steps
npoint=  100;  % number of vertices

% * create element decomposition for FE method
x                = linspace(0,L,npoint)';
nelement         = npoint-1;
nd(1:nelement,1) = 1:npoint-1;
nd(1:nelement,2) = 2:npoint;
local_mass_p1    = [1/3 1/6;1/6 1/3];

% * create & remember initial data
h  = L/2-abs(L/2-x); 
hi = h;xi = x;
t  = 0;dt = T/nt;

% * exact stationary solution
a  = sqrt(sqrt(18)/16);
xs = linspace(-a,a,100);
hs = a/sqrt(2)*(1-(xs/a).^2);xs=xs+1/2;
for it=1:nt       
    % * plot numerical solution with initial data and 
    % stationary solution for provided data
    if mod(it,10)==1
        plot(x,h);
        xlabel('x'); 
        ylabel('h');
        hold on
        drawnow
    end
    % * construct system matrices
    build_FE_matrices % script: matrices A,S,Ms,Dx for FEM
    build_ALE_matrix  % script: matrix for ALE decomposition
    
    % * FE problem: build right-hand-side rhs & solve
    rhs=[zeros(npoint,1);S*h+Ms*(c1*h+c2*x)];
    rhs(ndof+1)=rhs(ndof+1)+(SL+(dh(  1)^2)/2)/dh(  1);
    rhs(2*ndof)=rhs(2*ndof)-(SR+(dh(end)^2)/2)/dh(end);
    hdot = A\rhs; % solve for u=(hdot,pi)^t
    
    % * perform ALE decomposition & update solution
    udec  = P\hdot(1:ndof); % decompose hdot
    h(2:ndof-1)=h(2:ndof-1)+dt*udec(2:ndof-1);% update h
    x          =x          +dt*I*udec;        % update x

end
hold off