function [CDF CDFa] = RosenblattCDF(z,D,varargin)
% Returns the values of the CDF of the Rosenblatt distribution for values
% in the vector Z.  

if isempty(varargin)
    M = 50;
    N = 5;
    func = 'pdf';
else
    M = varargin{1};
    N = varargin{2};
    func = varargin{3};
end

sizeZ = size(z);
z = z(:);
CDF = zeros(size(z));
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

% CDF of sum of chi-squared distributions
F = @(x) real( (x > -mu+.01) .* ...
    nninfdiv(x + mu,{'chi-squared',1,Evals(:)},func));

% Compute CDF by approximating convolutions
for i = 1:length(z)
    % Approximate CDF until convergence
    CDFaOld = -inf;
    for n=N
        %Fconv = @(y) F(z(i) - y) .* phi(y/sigM) .* p{n}(y/sigM) / sigM;
        Fconv = @(y) f_convolution(y, F, z(i), phi, sigM, p{n});
        CDFa(n) = quadgk(Fconv,-inf,0) + quadgk(Fconv,0,mu);
        if abs(CDFa(n) - CDFaOld) < 1e-6 
            break
        end
        CDFaOld = CDFa(n);
    end
    CDF(i) = CDFa(n);
end

CDF = reshape(CDF,sizeZ);

end

function f = f_convolution(y, F, z, phi, sigM, p)
    f = F(z - y) .* phi(y/sigM) .* p(y/sigM) / sigM;
    f( isnan(f) ) = 0;
end



function f = hurwitz(s,M)

f = zeta(s) - sum((1:(M-1)).^(-s));

end







