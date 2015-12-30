function F = hypgeof(a,b,c,z)
% Sums the hypergeometric series 
% sum( (a)_n (b)_n / (n! (c)_n) z^n ) 
% Only valid for |z| < 1, and 
% for best results, assumes c > a+b and c > a*b


F = zeros(size(z));
mask = (z<=.5);

if any(mask)    
    
    M = 60; % number of terms in series to sum
    Q = zeros(1,M); 
    z1 = z(mask);
    t0 = [a b c 1];
    t = t0;
    Q(1) = 1;
    Q(2) = t(1) * t(2) / (t(3) * t(4));
    % Compute coefficients
    for i=1:M-2
        t = t .* (t0 + i);
        Q(i+2) = t(1) * t(2) / (t(3) * t(4));
    end
    Q = Q(end:-1:1);
    F(mask) = polyval(Q,z1);

end


if any(~mask)
    z1 = z(~mask);
    F(~mask) = gamma(c)*gamma(c-a-b)/(gamma(c-a)*gamma(c-b)) *...
                                            hypgeof(a,b,a+b+1-c,1-z1) ...
        + gamma(c) *gamma(a+b-c)/(gamma(a)*gamma(b)) *...
                       (1-z1).^(c - a - b) .* hypgeof(c-a,c-b,1+c-a-b,1-z1);
end
                                    



end




