function x = meshgen(a,b,N)
    x = linspace(a,b,N);
    y = x;
    y(x<=1/2) = x(x<=1/2).^4 * .5^(-3);
    y(x>1/2) = 1 - (1 - x(x>1/2)).^4 * .5^(-3);
    x = y;
end