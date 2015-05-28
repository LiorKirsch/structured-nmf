function [U,V] = init_with_svd(X,k, check_for_most_zeros)

if ~exist('check_for_most_zeros','var')
    check_for_most_zeros = true;
end

[U,S,V] = svds(X,k) ;

% U2 = max(U,0);
% V2 = max(V,0);

% This part flips the the vectors which are mostly zeros.
if check_for_most_zeros
   for i = 1:k
      if sum(V(:,i)  >0)+ sum(U(:,i) >0 ) < sum(V(:,i) < 0)+ sum(U(:,i) <0 )
          V(:,i) = -1 * V(:,i);
          U(:,i) = -1 * U(:,i);
      end
   end
end

U = max(U,0);
V = max(V,0);


fprintf('distance from svd init %g\n', norm( X - U*S*V' , 'fro') );
% fprintf('distance from svd init %g\n', norm( X - U2*S*V2' , 'fro') );
V = V*S;
V = V';

end