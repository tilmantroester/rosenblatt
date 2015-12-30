function CDF = RosenblattCDF(z,D,M,N)
% Returns the values of the CDF of the Rosenblatt distribution parametrized by 
% D for values in the vector z.
%
% Inputs:
%      z:    Value or array of values to compute CDF.
%      D:    Parameter D in (0,1/2) for Rosenblatt distribtion  
%      M:    (Optional) Number of chi-squared distributions to include in approximation
%      N:    (Optional) Number of terms to use in Edgeworth expansion to approximate infinite sum of chi-squares
%
% Outputs:  
%      CDF:   value of cdf
%
%  Note:  M,N take default values of 50 and 5, respectively.  To obtain more digits of accuracy 
%         it might be necessary to increase these values.  
%
%  Examples:  
%       cdf  =  RosenblattCDF( [0 1 2] , 0.3 );
%       cdf2 =  RosenblattCDF( [0 1 2] , 0.3, 100, 6 ); % Increases M,N
%
% Reference:
%   Veillette, Mark S., and Murad S. Taqqu. "Properties and numerical evaluation of the Rosenblatt distribution." Bernoulli 19.3 (2013): 982-1005.
%
%
if nargin == 2
    M = 50;
    N = 5;
end
CDF = rosenblatt_dist( z , D, 'cdf', M , N );

end


