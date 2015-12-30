function mu = c2m(kap)
% Transforms the first n cumulants into the first n moments.
% KAP is a vector of the first n cumulants.
% Returns a column vector MU of the first n moments.

% See P. J. Smith, "A recursive formulation of the old prolem of obtaining
% moments from cumulants and vice versa"  American Statistician, Vol. 49
% no. 2 1995 pp 217-218

kap = kap(:)';
N = length(kap);
mu = zeros(1,N);
mu(1) = kap(1);
pTriangle = zeros(N+1);
pTriangle(1,2) = 1;
for n = 2:N
    pTriangle(n,2:end) = pTriangle(n-1,1:end-1) + pTriangle(n-1,2:end);
    mu(n) = sum( pTriangle(n,2:n+1) .* kap(n:-1:1) .* [1 mu(1:n-1)] );
end
mu = mu(:);


end
