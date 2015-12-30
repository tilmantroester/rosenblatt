function y = rosenblatt_dist( z, D, func, M, N )
% Computing rosenblatt PDF/CDF of Rosenblatt distribution with parameter D.
% Inputs:
%      z:    Value or array of values to compute CDF/PDF.
%      D:    Parameter D in (0,1/2) for Rosenblatt distribtion 
%      func: String that is either 'pdf' or 'cdf'.  
%      M:    (Optional) Number of chi-squared distributions to include in approximation
%      N:    (Optional) Number of terms to use in Edgeworth expansion to approximate infinite sum of chi-squares
%
% Outputs:  
%      y:   value of pdf/cdf
%
%  Note:  M,N take default values of 50 and 5, respectively.  To obtain more digits of accuracy 
%         it might be necessary to increase these values.  
%
%  Examples:  
%       pdf  =  rosenblatt_dist( 1, 0.3, 'pdf' );
%       cdf  =  rosenblatt_dist( [0 1 2] , 0.3, 'cdf' );
%       cdf2 =  rosenblatt_dist( [0 1 2] , 0.3, 'cdf', 100, 6 );
%
% Reference:
%   Veillette, Mark S., and Murad S. Taqqu. "Properties and numerical evaluation of the Rosenblatt distribution." Bernoulli 19.3 (2013): 982-1005.
%

if nargin == 3
    M = 50;
    N = 5;
end

sizeZ = size(z);
z = z(:);
y = zeros(size(z));
conv = true(size(z));

% Find eigenvalues and correction terms
sig = sqrt(.5*(1-2*D)*(1-D));
CD = 2/(pi^(1-D)) * sig * gamma(1-D) * sin(pi*D/2);

% Find eigenvalues
Evals = RosenblattEigs(D,M,1500) * sig;
sigM = sqrt(1 - 2 * sum(Evals.^2));
%kap = RosenblattC(D,N,'cumulants');
kapM = zeros(N,1);
kapM(2) = 1;

for i = 3:N
        %kapM(i) = kap(i-1) - 2^(i - 1) * factorial(i-1) * sum(Evals.^i);
        kapM(i) = 2^(i-1)*factorial(i-1)*CD^i*hurwitz(i*(1-D),M+1);
        kapM(i) = kapM(i)/sigM^i;
end

% Standard normal pdf
c = 1/sqrt(2 *pi);
phi = @(x) c * exp(-.5*x.^2);

% Hermite polynomials
H = cell(3*(N-2),1);
for k = 1:3*(N-2)
    v = hermiteh(k);
    H{k} = @(x) polyval(v,x);
end

% Zeta + 1
zeta1 = @(kvec) sum(  (3:(length(kvec)+2)) .* kvec ); 

% Edgewood corrections for PDF of YM
p = cell(N,1);
p{2} = @(x) 1;
for n=3:N
    m = 3:n;
    p{n} = @(x) p{n-1}(x);
    etan = eta(n);
    for j=1:size(etan,1)
        k = etan(j,:);
        p{n} = @(x)  p{n}(x) +...
            prod(1./factorial(k) .* (kapM(m)'./factorial(m)).^k).*H{zeta1(k)}(x);
    end
end

% Shift
mu = sum(Evals);

% CDF/PDF of sum of chi-squared distributions
F = @(x) real( (x > -mu+.01) .* ...
    nninfdiv(x + mu,{'chi-squared',1,Evals(:)},func));

% Compute CDF by approximating convolutions
for i = 1:length(z)
    % Approximate CDF until convergence
    yaOld = -inf;
    for n=N
        %Fconv = @(y) F(z(i) - y) .* phi(y/sigM) .* p{n}(y/sigM) / sigM;
        Fconv = @(y) f_convolution(y, F, z(i), phi, sigM, p{n});
        ya(n) = quadgk(Fconv,-inf,0) + quadgk(Fconv,0,mu);
        if abs(ya(n) - yaOld) < 1e-6 
            break
        end
        yaOld = ya(n);
    end
    y(i) = ya(n);
end

y = reshape(y,sizeZ);

end

% Utility functions
function f = f_convolution(y, F, z, phi, sigM, p)
    f = F(z - y) .* phi(y/sigM) .* p(y/sigM) / sigM;
    f( isnan(f) ) = 0;
end

function f = hurwitz(s,M)
    f = zeta(s) - sum((1:(M-1)).^(-s));
end


