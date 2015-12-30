function PDF = RosenblattPDF(z,D,M,N)
% Returns the values of the PDF of the Rosenblatt distribution parametrized by 
% D for values in the vector z.
%
% Inputs:
%      z:    Value or array of values to compute PDF.
%      D:    Parameter D in (0,1/2) for Rosenblatt distribtion  
%      M:    (Optional) Number of chi-squared distributions to include in approximation
%      N:    (Optional) Number of terms to use in Edgeworth expansion to approximate infinite sum of chi-squares
%
% Outputs:  
%      PDF:   value of pdf
%
%  Note:  M,N take default values of 50 and 5, respectively.  To obtain more digits of accuracy 
%         it might be necessary to increase these values.  
%
%  Examples:  
%       pdf  =  RosenblattPDF( [0 1 2] , 0.3 );
%       pdf2 =  RosenblattPDF( [0 1 2] , 0.3, 100, 6 ); % Increases M,N
%
% Reference:
%   Veillette, Mark S., and Murad S. Taqqu. "Properties and numerical evaluation of the Rosenblatt distribution." Bernoulli 19.3 (2013): 982-1005.
%
%
if nargin == 2
    M = 50;
    N = 5;
end
PDF = rosenblatt_dist( z , D, 'pdf', M , N );

end


