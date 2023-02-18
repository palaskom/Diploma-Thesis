function [MI,g,p] = Descrete_Mutual_Information(Nmp,P,r,sigma)

% Nmp: number of mass points
% P: average power
% r = A/P
%     if (10*log10(P/sigma)<10)
%         N = 20;
% %         (10*log10(P/sigma)>10 || 10*log10(P/sigma)<5)
% %         N = 20;
%     else 
%         N = 30;
%     end


%     if (Nmp>1)
%         mean = (0:Nmp-1)*r*P/(Nmp-1);
%     else
%         mean = 0;
%     end
    if (Nmp<2)
        Nmp = 2;
    end
    
    if (Nmp==2)
        if (r<=2)
            p = [0.5 0.5];
        else
            p = [0.75 0.25];
        end
    else
        p = DiscretePDF(Nmp,r);
    end
    
    mean = (0:Nmp-1)*r*P/(Nmp-1);
    
%     x = mean(1)-4*sigma:1/N:mean(end)+4*sigma;
    if (10*log10(P/sigma)<=10)
        N = 500;
    else
        N = 300;
    end
    x = linspace(mean(1)-4*sigma,mean(end)+4*sigma,N);
    dx = x(2)-x(1);
    
%     Dx = mean(end)+4*sigma - ( mean(1)-4*sigma );
%     x = mean(1)-4*sigma:Dx/N:mean(end)+4*sigma;
%     dx = x(2)-x(1);
    

    f = 0;
    for j = 1:length(p)
        f = f + p(j)*(1/sqrt(2*pi*sigma^2))*exp(-(x-mean(j)).^2/(2*sigma^2));
    end

    g = f.*log2(f);
    MI = -sum(g)*dx - (1/2)*log2(2*pi*exp(1)*sigma^2);
        
end

