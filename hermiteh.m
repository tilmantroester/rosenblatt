function h = hermiteh(n)
%Returns a vector corresponding to the nth hermite polynomial
if n==0
    %h = @(x) 1;
    h = 1;
else 
    syms x y;
    y = x;
    for i = 2:n
        y = simplify(-diff( y * exp(-x.^2/2) ) * exp(x.^2/2)); 
    end
    %h = @(x) subs(y);
    h = sym2poly(y);
    
end