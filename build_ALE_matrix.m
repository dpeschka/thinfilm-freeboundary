% build matrix for decomposition of Eulerian time-derivative to ALE
% time-derivative and transport velocity
P  = sparse(ndof,ndof);
I  = sparse(ndof,ndof);
X  = sparse(ndof,ndof);

dh = M\(Dx*h);
xi = (x-x(1))/(x(end)-x(1));

for i=1:ndof
    P(i,1)     =( 1-xi(i)) * dh(i);
    P(i,npoint)=(   xi(i)) * dh(i);
    X(i,1)     =( 1-xi(i)) ;
    X(i,npoint)=(   xi(i)) ;
end
for i=2:npoint-1
    I(i,i)=1;
end