function [Rmin,v] = CapacityLowerBound(r,P,sigma)

    if (r<=2)
        v = r;
    else
        g = @(x) ( 1/x+1/(1-exp(x))-1/r );
        mu = fzero(g,1);
        v = r*(1-exp(-mu))/(mu*exp(-mu/r));
    end
    
%     Rmin = (1/2)*log(1+(v*P)^2/(2*pi*exp(1)*sigma^2));
    Rmin = (1/2)*log2(1+(v*P)^2/(2*pi*exp(1)*sigma^2));

end