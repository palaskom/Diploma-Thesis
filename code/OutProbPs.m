function [Po,fPS] = OutProbPs(P, r, d, rth, qth, sadc, srec, a, b, model)

% r = A/P
% P: average power constraint
% d: distance from user
% rth: rate threshold
% qth: energy threshold
% sadc: ADC noise
% srec: REC noise
% h~gamma(nT,sh^2), a=nT, b=sh^2
% model = 1 for the linear EH model and 2 for the nonlinear EH model

    % constants
    at = 0.5; % aperture of trasmit antenna [m]
    ar = 0.01; % aperture of receive antenna [m]
    fc = 2.45e9; % operating frequency
    L = 1-exp(-at*ar/d^2/(3e8/fc)^2); % path loss factor
    [~,v] = CapacityLowerBound(r,P,sadc);
    
    if (r<=2)
        E = r*P/2;
    else
        E = P;
    end
    
    % constant related to the harvested energy
    if (model==1)
        z = 0.78;
        c1 = qth/(z*L*E);
    elseif (model==2)
        Ps = 0.024; A = 150; B = 0.014;
        c1 = (1/A/L/E)*log( (Ps+qth*exp(A*B))/(Ps-qth) );
    end
    
    % constant related to the rate
    c2 = (1/v/L/P)*sqrt(2*pi*exp(1)*(2^(2*rth)-1));
    
    % polynomial coefficients
    pc(1) = 1;  pc(2) = -2;
    pc(3) = 1 + (sadc/srec)^2 - (c1/c2/srec)^2;
    pc(4) = 2*(c1/c2/srec)^2;
    pc(5) = -(c1/c2/srec)^2;
    
    % optimal value of PS factor
    root = roots(pc); % Get all roots
    fPS = root(imag(root)==0 & real(root)>0 & real(root)<=1); % Save only the real roots @ (0,1]
    
    % outage probability
    k1 = c1/fPS; 
    k2 = c2*sqrt(srec^2 + sadc^2/(1-fPS)^2);
    Po = gammainc(max(k1,k2)/b,a);
    
end

