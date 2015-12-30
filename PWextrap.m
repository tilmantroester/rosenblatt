function f = PWextrap(k,P,method)
% Performs Post-Widder extrapolation to approximate convergence of a finite sequence.
%  Inputs:  
%      k:  Sequence of indices in sequence 
%         (for 'polynomial', this must be in the form k_j = k_0 * j for some k_0)
%      P:  Terms in sequence whose limit is being approximated
%      method:   Extrapolation method:  'polynomial' or 'rational'. 
%  Output:  
%      f :  Approximation of limit of sequence.
%
%
    switch lower(method)
        case 'polynomial'
            % Extrapolates limit using Lagrange polynomial interpolation
            N = size(P,3);
            j = length(k);
            f = zeros(1,1,N);
            c = zeros(1,j,N);
            % Assumes k_j = k_0 * j !!!!!
            for r=1:j
                c(1,r,:) = repmat( ...
                    1./( prod(1 - ( (1:(r-1))/ r ) ) * ...
                    prod(1 - ((r+1):j)/r ) )...
                            ,[1,1,N]);
                f = f + P(1,r,:).*c(1,r,:);
            end
                 
        case 'rational'
            % Extrapolates limit using rational function interpolation
            r = size(P,2);
            d = size(P,3);
            T = zeros(r+1,r,d);
            T(2,:,:) = P(1,end:-1:1,:);
            k = k(end:-1:1);
            for i=3:(r+1)
                for j=1:(r-i+2)
                    T(i,j,:) = T(i-1,j+1,:) + (T(i-1,j+1,:)-T(i-1,j,:))...
                               ./ ( ( k(i+j-2)/k(j) )  ...
                               .* ( T(i-1,j,:) - T(i-2,j+1,:) )...
                               ./ ( T(i-1,j+1,:) - T(i-2,j+1,:) ) - 1 );
                end
            end

            f = T(end,1,:);
     end


end