% build matrix for decomposition of Eulerian time-derivative to ALE
% time-derivative and transport velocity
P  = sparse(ndof,ndof);
I  = sparse(ndof,ndof);
X  = sparse(ndof,ndof);

dh = M\(Dx*h);
xi = (x-x(1))/(x(end)-x(1));

iip = [1:ndof 1:ndof];
jjp = [1+zeros(1,ndof) npoint+zeros(1,ndof)];
iii = 2:npoint-1;

P = sparse(iip,jjp,[(1-xi).*dh xi.*dh],ndof,ndof);
X = sparse(iip,jjp,[(1-xi) xi],ndof,ndof);
I = sparse(iii,iii,ones(1,npoint-2),ndof,ndof);

% for i=1:ndof
%     P(i,1)     =( 1-xi(i)) * dh(i);
%     P(i,npoint)=(   xi(i)) * dh(i);
%     X(i,1)     =( 1-xi(i)) ;
%     X(i,npoint)=(   xi(i)) ;
% end
% for i=2:npoint-1
%     I(i,i)=1;
% end