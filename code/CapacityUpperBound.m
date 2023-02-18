function Rmax = CapacityUpperBound(r,P,sigma)

    A = r*P;
    a = 1/r;
    delta = sigma*log(1+A/sigma);    

    if (r<=2)
        Rmax = (1-2*qfunc((delta+A/2)/sigma)).*log( (A+2*delta)./(sigma*sqrt(2*pi)*(1-2*qfunc(delta/sigma))) )...
               -1/2 + qfunc(delta/sigma) + (delta/sigma/sqrt(2*pi)).*exp(-delta.^2/2/sigma^2);
    else
        g = @(x) ( 1/x+1/(1-exp(x))-a );
        m = fzero(g,1);
        mu = m*(1-exp(-a*delta^2/2/sigma^2));
        
        Rmax = ( 1-qfunc((delta+a*A)/sigma)-qfunc((delta+(1-a)*A)/sigma) )...
           *log( (A/sigma) * ( exp(mu*delta/A)-exp(-mu*(1+delta/A)) )/( sqrt(2*pi)*mu*(1-2*qfunc(delta/sigma)) ) )...
           -1/2 + qfunc(delta/sigma) + ( delta/(sigma*sqrt(2*pi)) )*exp(-delta^2/2/sigma^2)...
           +(mu*sigma/A/sqrt(2*pi))*( exp(-delta^2/2/sigma^2) - exp(-(A+delta)^2/2/sigma^2) )...
           +mu*a*( 1-2*qfunc((delta+A/2)/sigma) );
    end
        
end

