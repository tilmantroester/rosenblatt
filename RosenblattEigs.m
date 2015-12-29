function [Lam phis T] = RosenblattEigs(D,M,varargin)
% Returns approximations of the frist M eigenvalues of the Rosenblatt
% distribution with parameter D.

if M > 5000
    error('Make M smaller!')
end

% Define discrete integral operator
% a = 0;
% b = 1;
% N = 5000;
% x = linspace(a,b,N);
% row1 = zeros(N,1);
% row1(1) = N^D/(1-D);
% row1(2:end) = abs(x(1) - x(2:end)).^(-D);
% T = toeplitz(row1)/N;
% T(:,1) = .5 * T(:,1);
% T(:,end) = .5 * T(:,end);
% Lam = eigs(T,M);


a = 0;
b = 1;
if isempty(varargin)
    N = 600;
else
    N= varargin{1};
end
x = meshgen(a,b,N);
%x = linspace(a,b,N);
e = x(2:end)-x(1:end-1);
T = zeros(N);
c = zeros(1,N-1);
d = zeros(1,N-1);
for i=1:N
    jlower = 1:(i-1);
    jupper = i:(N-1);
    xlower = x(jlower);
    xlowerp1 = x(jlower + 1);
    xupper = x(jupper);
    xupperp1 = x(jupper + 1);
    if any(jlower)  
        c(jlower)= (x(i) - xlower).^(1-D) / (1-D) -  ...
                    ((x(i) - xlower).^(2-D) - (x(i) - xlowerp1).^(2-D))...
                    ./ ( (1-D) * (2-D) * e(jlower) );
        
        d(jlower)=( (x(i) - xlower).^(2-D) - (x(i) - xlowerp1).^(2-D)) ...
                  ./( (1-D)*(2-D)*e(jlower) ) ...
                  - (x(i) - xlowerp1).^(1-D) / (1-D);
    
    end
    
    if any(jupper)
        c(jupper)= ((xupperp1 - x(i)).^(2-D) - (xupper - x(i)).^(2-D))...
                   ./( (1-D)*(2-D)*e(jupper) ) ...
                   - (xupper - x(i)).^(1-D) / (1-D);
        
        d(jupper)= (xupperp1 - x(i)).^(1-D) / (1-D) ...
                   -((xupperp1 - x(i)).^(2-D) - (xupper - x(i)).^(2-D))...
                   ./( (1-D)*(2-D) *e(jupper) );
    
    end
    T(i,1) = c(1);
    T(i,end) = d(end);
    T(i,2:end-1) = c(2:end) + d(1:end-1);
    
end


if nargout <=1
    Lam = eigs(T,M);
elseif nargout >1
    [phis Lam] = eigs(T,M);
    Lam = Lam(Lam~=0);
end


end









