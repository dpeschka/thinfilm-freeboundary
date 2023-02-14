xii=([-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))]+1)/2;
wii=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]/2;
for k=1:nelement
%     dphi = [-1 1]/edet(k);         % local gradient   
%     % build local matrices    
%     sloc = (dphi'*dphi)  * edet(k);% stiffness matrix
%     mloc = local_mass_p1 * edet(k);% mass matrix
%     cloc = [dphi;dphi]/2 * edet(k);% derivative matrix
%     % generate indices for sparse matrix
     ii(k,:) = [nd(k,1) nd(k,2) nd(k,1) nd(k,2)];
     jj(k,:) = [nd(k,1) nd(k,1) nd(k,2) nd(k,2)];    
%     % generate values for sparse matrix
%     S_ (k,:) =         sloc(:);
%     Ms_(k,:) =         mloc(:);
%     Sw_(k,:) = mobi(k)*sloc(:);
%     Dx_(k,:) =         cloc(:);
end