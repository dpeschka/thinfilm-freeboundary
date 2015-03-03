% build matrix for decomposition of Eulerian time-derivative to ALE
% time-derivative and transport velocity
P  = sparse(ndof,ndof);
I  = sparse(ndof,ndof);
dh = Ms\(Dx*h);
xi = (x-x(1))/(x(end)-x(1));
for i=1:ndof
    P(i,1)     =-( 1-xi(i)) * dh(i);
    P(i,npoint)=-(   xi(i)) * dh(i);
    I(i,1)     =+( 1-xi(i)) ;
    I(i,npoint)=+(   xi(i)) ;
end
for i=2:npoint-1
    P(i,i)=1;
end