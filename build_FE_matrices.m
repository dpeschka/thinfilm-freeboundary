% build matrices for the finite element method
ndof = npoint;                % number of degrees of freedoms
edet = x(nd(:,2))-x(nd(:,1)); % determinant of transformation
ii  = zeros(nelement,4); % integer array of indices
jj  = zeros(nelement,4); % integer array of indices
S_  = zeros(nelement,4); % real  array of matrix values
Ms_ = zeros(nelement,4); % real  array of matrix values
Sw_ = zeros(nelement,4); % real  array of matrix values
Dx_ = zeros(nelement,4); % real  array of matrix values

% explicit integration for mobility with alpha=2
mobi = zeros(nelement,1);
for i=1:2
    for j=1:2
        mobi = mobi + local_mass_p1(i,j)*h(nd(:,i)).*h(nd(:,j));
    end
end

% build global matrices from local matrices 
for k=1:nelement    
    dphi = [-1 1]/edet(k);         % local gradient   
    % build local matrices    
    sloc = (dphi'*dphi)  * edet(k);% stiffness matrix
    mloc = local_mass_p1 * edet(k);% mass matrix
    cloc = [dphi;dphi]/2 * edet(k);% derivative matrix
    % generate indices for sparse matrix
    ii(k,:) = [nd(k,1) nd(k,2) nd(k,1) nd(k,2)];
    jj(k,:) = [nd(k,1) nd(k,1) nd(k,2) nd(k,2)];    
    % generate values for sparse matrix
    S_ (k,:) =         sloc(:);
    Ms_(k,:) =         mloc(:);
    Sw_(k,:) = mobi(k)*sloc(:);
    Dx_(k,:) =         cloc(:);
end
% generate sparse matrices, e.g. S(ii(k),jj(k))=S_(k)
S  = sparse(ii(:), jj(:), S_(:)  );
Sw = sparse(ii(:), jj(:), Sw_(:) );
Ms = sparse(ii(:), jj(:), Ms_(:) );
Dx = sparse(ii(:), jj(:), Dx_(:) );
Z  = sparse(ndof,ndof);
% build 2x2 block matrices from ndof x ndof matrices
A = [Ms Sw;-dt*S Ms];
M = [Ms  Z;Z  Ms];