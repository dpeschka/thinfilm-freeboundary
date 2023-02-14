% build matrices for the finite element method
%integration_weights;

ndof = npoint;                % number of degrees of freedoms
edet = x(nd(:,2))-x(nd(:,1)); % determinant of transformation
%ii  = zeros(nelement,4); % integer array of indices
%jj  = zeros(nelement,4); % integer array of indices
%S_  = zeros(nelement,4); % real  array of matrix values
%Ms_ = zeros(nelement,4); % real  array of matrix values
%Sw_ = zeros(nelement,4); % real  array of matrix values
%Dx_ = zeros(nelement,4); % real  array of matrix values

% 5-point Gauss integration
mobi = zeros(nelement,1);
for i=1:length(xii)
    
    hval = (1-xii(i))*h(nd(:,1)) + xii(i) * h(nd(:,2));    
    mobi = mobi + wii(i) * abs(hval).^mobexp;
    
end


% % explicit integration for mobility with alpha=2
% mobi = delta*mean(h(nd),2);%zeros(nelement,1);
% mobi = mobi + alpha;
% for i=1:2
%     for j=1:2
%         mobi = mobi + beta*local_mass_p1(i,j)*h(nd(:,i)).*h(nd(:,j));
%     end
% end
% 
% for i=1:2
%     for j=1:2
%         for k=1:2
%             mobi = mobi + mu/3*local_cube_p1(i,j,k)*h(nd(:,i)).*h(nd(:,j)).*h(nd(:,k));
%         end
%     end
% end

dphiref = [-1 1];
% build global matrices from local matrices
% for k=1:nelement
%     dphi = [-1 1]/edet(k);         % local gradient   
%     % build local matrices    
%     sloc = (dphi'*dphi)  * edet(k);% stiffness matrix
%     mloc = local_mass_p1 * edet(k);% mass matrix
%     cloc = [dphi;dphi]/2 * edet(k);% derivative matrix
%     % generate indices for sparse matrix
%     ii(k,:) = [nd(k,1) nd(k,2) nd(k,1) nd(k,2)];
%     jj(k,:) = [nd(k,1) nd(k,1) nd(k,2) nd(k,2)];    
%     % generate values for sparse matrix
%     S_ (k,:) =         sloc(:);
%     Ms_(k,:) =         mloc(:);
%     Sw_(k,:) = mobi(k)*sloc(:);
%     Dx_(k,:) =         cloc(:);
% end
% vectorized version
S_   = zeros(nelement,2,2);
Ms_  = zeros(nelement,2,2);
Dx_  = zeros(nelement,2,2);
Sw_  = zeros(nelement,2,2);

for i=1:2
    for j=1:2
        S_(:,i,j)  = S_(:,i,j)  + dphiref(i)*dphiref(j)./edet;
        Sw_(:,i,j) = Sw_(:,i,j) + dphiref(i)*dphiref(j).*mobi./edet;       
        Ms_(:,i,j) = Ms_(:,i,j) + local_mass_p1(i,j) .* edet;
        Dx_(:,i,j) = Dx_(:,i,j) + dphiref(j)/2;
    end
end

% generate sparse matrices, e.g. S(ii(k),jj(k))=S_(k)
S  = sparse(ii(:), jj(:), S_(:)  );
Sw = sparse(ii(:), jj(:), Sw_(:) );
M  = sparse(ii(:), jj(:), Ms_(:) );
Dx = sparse(ii(:), jj(:), Dx_(:) );
Z  = sparse(ndof,ndof);
% build 2x2 block matrices from ndof x ndof matrices
A = [M  Sw;-dt*S M];