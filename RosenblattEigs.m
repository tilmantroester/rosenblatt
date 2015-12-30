function [Lam phis T] = RosenblattEigs(D,M,N)
% Returns approximations of the first M eigenvalues of the Rosenblatt
% distribution with parameter D.
%  Inputs:
%      D:  Parameter in Rosenblatt distribution in (0,1/2)
%      M:  Number of eigenvalues to return (must be less than N)
%      N:  (Optional) Number of gridpoitns used to approximate discrete integral operator
%          Defaults to 600.
%  Note:  N should be much bigger than M.
%
%  Outputs:
%      Lam:  Vector of M eigenvalues
%     phis:  Matrix of eigenvectors
%     T   :  Discrete matrix operator using to compute eigenvalues
%
%  Example:
%      eig_vals = RosenblattEigs(0.3)
%

if nargin == 2
    N = 600;
end

if M > 5000
    error('You should probably make M smaller.  Comment this error out if you''re feeling lucky..')
end

if M > N
    error('M must be less or equal to N');
end

a = 0;
b = 1;
x = meshgen(a,b,N);
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









