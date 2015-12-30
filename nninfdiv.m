function varargout = nninfdiv(x,dist,varargin)
% F = nninfdiv(X,DIST) returns the values of the PDF of the non-negative
% infinitely divisible distribution described in DIST at the values 
% in X.  DIST must cel a cell array which has the following format
%   DIST = {NAME,PARAM1,PARAM2,...}
% NAME is a string which is the name of the dsitribution.
% Possible options for DIST are
%
%   DIST = {'alpha stable',[alpha1 alpha2... alphaN],[c1 c2...cN]}
%       Gives sum of N alpha stable distributions with alpha's
%       alpha1...alphaN weighted by non-negative constants c1...cN.  
%
%   DIST = {'uniform mix'}
%       Uniform mixture of alpha stable multiplied by C>0 
%
%   DIST = {'chi-squared',K,C}
%       Sums of Chi-squared random with K degrees of freedom.
%       K is a scalar or vector.  C is a scalar or vector.
%       gives pdf/cdf of sum C(1) X(1) + C(2) X(2) + ... + C(N) X(N)
%       where X(1) has K(1) degress of freedom etc...
%       If K (C) is a scalar and C (K) is a vector, then replaces scalar
%       with vector same length of C whose entries are equal to K.
%
% Reference:
%     Veillette, Mark S., and Murad S. Taqqu. "A technique for computing the PDFs and CDFs of 
%       nonnegative infinitely divisible random variables." Journal of Applied Probability 48.1 (2011): 217-237.
%

% Set defaults
TOL = 1e-6;
func = 'pdf';
extrap = 'polynomial';
showcvg = false;
% Gather optional inputs
for i=1:length(varargin)
    if strcmpi(varargin{i},'cdf')
        func = 'cdf';
    elseif strcmp(varargin{i},'derivative')
        func = 'der';
        q = varargin{i+1};
        varargin{i+1} = false;
    elseif strcmpi(varargin{i},'rational')
        extrap = 'rational';
    elseif isnumeric(varargin{i})
        TOL = varargin{i};
    elseif strcmp(varargin{i},'converged')
        showcvg = true;
    end
    
end


% Deal with x
size_x = size(x); % so we can reshape later
x = x(:); % Make it a column vector
x3 = permute(x,[2 3 1]);
N = length(x);
% Keep track of error and x-values which haven't converged.
nocvg = true(N,1); 
abserr = inf(1,1,N);



Nks = 10; % Number of k's to use
k0 = 10;  % Use evenly spaced k's, k0, 2*k0, 3*k0 etc...

kvec = k0*(1:Nks); % Vector of k's
kmax = k0*Nks;

% Step 1. Initialize 
% Set up D matrices, 3rd dimension is for x values.
Dpsi = zeros(kmax,Nks,N);  % Holds derivatives of psi 
Dphi = zeros(kmax-1,Nks,N); % Holds derivatives of phi

% Vectors for approximations.
P = zeros(1,Nks,N);
f = zeros(1,Nks,N);

% Construct Pascal's triangle as we go for binomial coefficients
PTriangle = zeros(kmax+1);
PTriangle(1,2) = 1;
pindex = 1;


% Main loop
for j = 1:Nks    
    
    kj = kvec(j);    
    N = nnz(nocvg);
    
    % Steps 2
    Dpsi(1,j,nocvg) = IDPsi(kj,x(nocvg),dist);
    
    % Step 3.
    Dphi(1:(kj-1),j,nocvg) = IDPhik(kj,x3(nocvg),dist);
    
    for i=2:kj
        % If needed, compute next row in Pascal Triangle
        if pindex < i
            PTriangle(i,2:end) = PTriangle(i-1,1:end-1)+PTriangle(i-1,2:end);
            pindex = pindex + 1;
        end
        r = (1:i-1)';  
        binmat = PTriangle(i-1,2:i)';
        binmat = binmat(:,1,ones(N,1)); %(indexing trick to avoid repmat)
        Dpsi(i,j,nocvg) = sum( ...
                        binmat .* Dpsi(r,j,nocvg) .* Dphi(i-r,j,nocvg) ...
                           ,1) ;
    end
    
    
    switch func
        case 'pdf'
        % Step 4.    
        z = Dpsi(kj,j,nocvg);
        P(1,j,nocvg) = (-1)^(kj-1)/factorial(kj-1) * ...
                        sign(z) .* ...
                       exp(kj * log(kj./x3(nocvg)) + log(abs(z)) );                  
        case 'cdf' 
        % Step 4*.
        r = repmat((0:kj-1)',[1,1,N]);
        xx = repmat(x3(nocvg),[kj,1,1]);
        z = Dpsi(kj - (0:kj-1),j,nocvg);
        P(1,j,nocvg) = sum( (-1).^(kj +r -1) ./ factorial(kj - r  - 1) .*  ...
                   sign(z) .* exp((kj-r-1).*log(kj./xx) + log(abs(z))) ... 
                    ,1);
                
        case 'der'
        % "Other useful derivatives"
        r = repmat((0:min(q,kj-1))',[1,1,N]);
        xx = repmat(x3(nocvg),[size(r,1),1,1]);
        z = Dpsi(kj - (0:min(q,kj-1)),j,nocvg);  
        P(1,j,nocvg) = (-1)^(kj-1)*factorial(q) * ...
            sum(1./(factorial(kj-r-1).*factorial(r).*factorial(q-r)) ...
                   .* sign(z) .* ...
                   exp((kj+q-r) .* log(kj./xx) + log(abs(z)) ),1);
                          
    end              

    
    % Step 5.
    if j==1
        f(1,1,:) = P(1,1,:);
    else
        f(1,j,nocvg) = PWextrap(kvec(1:j),P(1,1:j,nocvg),extrap);
        f(1,j,~nocvg) = f(1,j-1,~nocvg);
        abserr(nocvg) = abs((f(1,j-1,nocvg) - f(1,j,nocvg))./f(1,1,nocvg));
        nocvg = abserr > TOL;
        if  all(~nocvg) || j == Nks
            p = reshape(f(1,j,:),size_x);
            break
        end   
    end
            
end

if strcmp(func,or('pdf','cdf'))
    p = max(0,p); % In case we extrapolated "too far".
end


varargout{1} = reshape(p,size_x);
if showcvg
    varargout{2} = ~nocvg(:);
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to compute Laplace transform psi(x) for various distributions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = IDPsi(k,x,dist)


    lambda = k./x;
    switch lower(dist{1})
    
        case {'alphastable','alpha stable'}
            
            alpha = dist{2};
            c = dist{3};
            logy = 0;
            for i=1:length(alpha)
                logy = logy - c(i)*lambda.^alpha(i) ;
            end
            y = exp(logy);          
            
        case {'uniform mix','uniformmix'}
            
            y = zeros(size(lambda));
            bad = abs(log(lambda)) < 1;
            % For lambda near 1, use Taylor series
            y(bad) = exp(-polyval(1./factorial(16:-1:1),log(lambda(bad))));
            % For lambda far from 1, use exact form
            y(~bad) = exp(-log(lambda(~bad)).^(-1) .* (lambda(~bad) - 1));
        
        
        case {'chi-squared','chisquared'}
                
            df = dist{2}; % degrees of freedom
            C = dist{3}; % Scales
            if isscalar(df) && isvector(C)
                df = df * ones(size(C));
            elseif isscalar(C) && isvector(df)
                C = C * ones(size(df));
            end
            
            y = 1;
            for i=1:length(C)    
                y = y .* (1 + 2 * C(i) * lambda).^(-df(i)/2);
            end
            
        
        case {'oupoisson','ou poisson'}    
            
            eta = dist{2};
            egamma = 0.5772156649015329;
            y = eta*(egamma + expint(lambda) + log(lambda));
            y = exp(-y);
            
            
            
        case {'ougamma','ou gamma'}
            eta = dist{2};
            kappa = dist{3};
            y = eta * kappa * dilog(1 + lambda);
            y = exp(-y);
            
            
            
        case {'custom'}
            phinder = dist{2};
            y = exp(-phinder(lambda,0));
            
            
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute derivatives of Laplace exponent 
%  for various distributions
%%%%%%%%%%%%%%%%%%%%%%%%%


function D = IDPhik(k,x,dist)


    switch lower(dist{1})
        
        case {'alphastable','alpha stable'}
        
            alpha = dist{2};
            c = dist{3};
            n = (1:(k-1))';
            Xmat = repmat(x,[k-1,1,1]);
            Nmat = repmat(n,[1,1,length(x)]);
            
            D = zeros(size(Xmat));
            for i=1:length(alpha)
                D = D + c(i) * (-1).^(Nmat+1) ...
                    .* (k./Xmat).^(alpha(i) - Nmat) ...
                    .* gamma(Nmat - alpha(i))/gamma(-alpha(i));
            end
  
        
        case {'uniformmix','uniform mix'}
            
            D = DphiUnifMix(k,x);
        
        case {'chi-squared','chisquared'}
        
            df = dist{2}; % degrees of freedom
            C = dist{3}; % Scales
            if isscalar(df) && isvector(C)
                df = df * ones(size(C));
            elseif isscalar(C) && isvector(df)
                C = C * ones(size(df));
            end
            
            
            n = (1:(k-1))';
            Xmat = repmat(x,[k-1,1,1]);
            Nmat = repmat(n,[1,1,length(x)]);

            D = zeros(size(Nmat));
            for i = 1:length(C)
            D = D + (-1).^Nmat .* (df(i)/2) ...
                .* exp(-Nmat.*log(1/(2*C(i)) + k./Xmat) + gammaln(Nmat));
            end
            
            
        case {'oupoisson','ou poisson'}  
            
            
            n = (1:(k-1))';
            Xmat = repmat(x,[k-1,1,1]);
            Nmat = repmat(n,[1,1,length(x)]);
            eta = dist{2};
            
            D = (-1).^(Nmat) * eta ./(k./Xmat).^(Nmat) ...
                .* gamma(Nmat).*gammainc(k./Xmat,Nmat);
            
            
        case {'ougamma','ou gamma'}
            
            n = 1:k-1;
            D = zeros(k-1,1,length(x));
            eta = dist{2};
            kappa = dist{3};
            lambda = k./x;
            
            for i=n
                D(i,1,:) = (-1).^i .* eta .* kappa .* ...
                    factorial(i-1)./lambda.^i .* (log(1 + lambda) ...
                    - polyval([1./(i-1:-1:1) 0],lambda./(1 + lambda)));
            end
            
            
    end


end



















