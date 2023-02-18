function p = DiscretePDF(Nmp,r)
% Nmp: number of mass points
% r = A/P
    
    K = Nmp - 1;
    
    if (r<=2)
        p = (1/Nmp)*ones(1,Nmp);
    else
        c = 1 - (K:-1:0)*r/K; % polynomial coefficients
        p = zeros(1,Nmp);

        root = roots(c); % Get all roots
        t = root(imag(root)==0 & real(root)>0 & real(root)<=1); % Save only the real roots @ (0,1]
        tsum = 0;
        for k = 0:K
            p(k+1) = t^k;
            tsum = tsum + p(k+1);
        end

        p = p/tsum;
    end
        
end

