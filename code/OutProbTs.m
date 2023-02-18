function [Po,fTS] = OutProbTs(P, r, d, rth, qth, sadc, srec, a, b, model)

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
    
    if (model==1)
        z = 0.78;
        fun = @(x) x^2*(2^(2*rth/(1-x))-1) - (v*qth*P)^2/(2*pi*exp(1)*(z*E)^2*(srec^2+sadc^2));
        fTS = fzero(fun,[0 0.99]);
        k1 = qth/(z*fTS*L*E);
    elseif (model==2)
        Ps = 0.024; A = 150; B = 0.014;
        fun = @(x) (2^(2*rth/(1-x))-1)/(log( (x*Ps-qth)/(x*Ps+qth*exp(A*B)) ))^2 - (v*P)^2/(2*pi*exp(1)*(A*E)^2*(srec^2+sadc^2));
        fTS = fzero(fun,[1e-4+qth/Ps 0.99]);
        k1 = (1/A/L/E)*log( (Ps*fTS+qth*exp(A*B))/(Ps*fTS-qth) );
    end
    
    k2 = (1/v/L/P)*sqrt( 2*pi*exp(1)*(srec^2+sadc^2)*(2^(2*rth/(1-fTS))-1) );
    Po = gammainc(max(k1,k2)/b,a);
    
end

