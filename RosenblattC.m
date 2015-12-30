function c = RosenblattC(D,K,type)
% Numerically approximates the first K cumulants of the Rosenblatt
% distribution with parameter D that is normalized to have variance 1.
%
% Inputs:
%    D:  Parameter of Rosenblatt distribution in (0,1/2)
%    K:  Number of cumulants to approximate
%    type:  Type of output.  Either 'cumulants', 'moments' or 'cs'.  
%    

if nargin == 2
    type = [];
 end

if D==0
    c = ones(max(K,3),1);
    sig = (2 * c(1) ).^(-.5);
elseif D<1/2
    c = zeros(max(K,3),1);
    c(1) = 1/((1-D)*(1 - 2*D));
    sig = (2 * c(1) ).^(-.5);
    c(2) = 2/((1-D)*(2-3*D)) * beta(1-D,1-D);

    % Define G1
    G1 = @(u) beta(1-D,1-D) * u.^(1-2*D) + ...
        (1-u).^(1-D) /(1-D) .* hypgeof(1,D,2-D,1-u);
    % Make accurate spline
    x = [0 meshgen(30,.005)];
    y1 = G1(x);
    G1sp = csapi(x,y1);

    % Compute c3 by integrating G1^2
    G1s = @(u) fnval(G1sp,u).^2;
    c(3) = (1/(1-D))*quadgk(G1s,0,1,'RelTol',1e-13);



    for k = 4:K    
        if mod(k,2)==0  % if even, we need to compute a new Gk
            y0 = Gk(x,G1sp,D);
            G0sp = csapi(x,y0);

            % Now compute cumulant
            F = @(u) fnval(G1sp,u).*fnval(G0sp,u);
            c(k) = (1/(1-D)) * quadgk(F,0,1);
        else         % if odd, we can reuse last Gk

            % Use Gk already computed
            F = @(u) fnval(G0sp,u).^2;
            c(k) = (1/(1-D)) * quadgk(F,0,1,'RelTol',1e-13);

            % reassign
            G1sp = G0sp; 
        end
    end

end

if ~isempty(type)
    switch type
        case 'cs'
            c = c;
        case 'cumulants'
            if D < 0.5
                Kv = (2:K+1)';
                c = 2.^(Kv-1) .* factorial(Kv - 1) .* sig.^(Kv) .* c;
            else
                c = zeros(max(K,3),1);
                c(1) = 1;
            end
        
        case 'moments'
            if D < 0.5
                Kv = (2:K+1)';
                c = 2.^(Kv-1) .* factorial(Kv - 1) .* sig.^(Kv) .* c;
            else
                c = zeros(max(K,3),1);
                c(1) = 1;
            end
            
            c = c2m([0 ;c]);
    end      
end





end


function m = meshgen(M,dx)
    m = exp(-M:dx:0);
    m = [m (1-exp(-(0:dx:M))/2)];
end


function y = Gk(u,pp,D)
    
    D1 = 1-D;
    
    Gexp1 = @(z) exp( D1 * z) .* fnval(pp,u.*(1-exp(z)));
    Gexp2 = @(z) exp( D1 * z) .* fnval(pp,u + exp(z)*(1-u));
    
    y = u.^(1-D) .* quadv(Gexp1,-30,0) ...
        + (1-u).^(1-D) .* quadv(Gexp2,-30,0);  
    
end










