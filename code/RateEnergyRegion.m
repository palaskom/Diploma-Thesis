%% channel
clear; clc;

simulations = 1; 

if (simulations==1)
    Niter = 1e3;
    nT = 4; % number of transmit antennas
    kapa = nT; % shape parameter 
    theta = 1; % scale parameter 
    h = gamrnd(kapa,theta,1,Niter);
else
    Niter = 1;
    h = 1;
end

%% constants

% path loss
at = 0.5; % aperture of trasmit antenna [m]
ar = 0.01; % aperture of receive antenna [m]
fc = 2.45e9; % operating frequency
d = 12; % distance in [m] between BS and the user
L = 1-exp(-at*ar/d^2/(3e8/fc)^2); % path loss factor

% EH model: 1 for linear EH and 2 for nonlinear EH
model = 2;
if (model==1)
    z = 0.78; % enenrgy conversion efficiency
else
    A = 150; B = 0.014; Ps = 0.024;
end

% SNR at the receiver
P = 2*L; % [W]
snr_dB = 10; % [dB]
snr = 10^(snr_dB/10);
s = P/snr;
r = 3; % r = A/P

% average transmit power
if (r>2)
    E = P;
else
    E = r*P/2;
end


%% RE Region - PS scheme

sf = 0:1e-2:1; % splitting factor

I_sum = zeros(1,length(sf));
Rmin_sum = zeros(1,length(sf));
Q_sum = zeros(1,length(sf));

for n = 1:length(h)    
    
    I = zeros(1,length(sf));
    Rmin = zeros(1,length(sf));
    Q = zeros(1,length(sf));
    
    Nmp = max(floor(h(n)*r*snr/2.5),2);
    for i = 1:length(sf)             
        sn = s*sqrt(1+1/(1-sf(i))^2); % standard deviation for PS 
        I(i) = Descrete_Mutual_Information(Nmp,h(n)*P,r,sn);
        Rmin(i) = CapacityLowerBound(r,h(n)*P,sn);
        if (model==1)
            Q(i) = z*sf(i)*h(n)*E;
        else
            Q(i)=(Ps/exp(A*B))*( (1+exp(A*B))/(1+exp(-A*(sf(i)*E*h(n)-B)))-1 );
        end

    end
    
    I_sum = I_sum + I;
    Rmin_sum = Rmin_sum + Rmin;
    Q_sum = Q_sum + Q;
              
end

I = I_sum/Niter;
Rmin = Rmin_sum/Niter;
Q = Q_sum/Niter;

I(end) = 0;
hold on, plot(Rmin,Q*1e3,'--g')
hold on, plot(I,Q*1e3,'g')

xlabel('Rate [bits/channel use]'), ylabel('Harvested Energy [mJ]')
legend('1x1 (R_{lb},Q_{NL})','1x1 (R_{ca},Q_{NL})','2x1 (R_{lb},Q_{NL})','2x1 (R_{ca},Q_{NL})',...
       '4x1 (R_{lb},Q_{NL})','4x1 (R_{ca},Q_{NL})')
grid on

%% RE Region - TS scheme

sf = 0:1e-2:1; % splitting factor

I_sum = zeros(1,length(sf));
Rmin_sum = zeros(1,length(sf));
Q_sum = zeros(1,length(sf));

for n = 1:length(h)    
    
    I = zeros(1,length(sf));
    Rmin = zeros(1,length(sf));
    Q = zeros(1,length(sf));
    
    Nmp = max(floor(h(n)*r*snr/2.5),2);
    for i = 1:length(sf)
              
        sn = s*sqrt(2); % standard deviation
        I(i) = (1-sf(i))*Descrete_Mutual_Information(Nmp,h(n)*P,r,sn);
        Rmin(i) = (1-sf(i))*CapacityLowerBound(r,h(n)*P,sn);
        if (model==1)
            Q(i) = z*sf(i)*h(n)*E;
        else
            Q(i)= sf(i)*(Ps/exp(A*B))*( (1+exp(A*B))/(1+exp(-A*(E*h(n)-B)))-1 );
        end
    end
    
    I_sum = I_sum + I;
    Rmin_sum = Rmin_sum + Rmin;
    Q_sum = Q_sum + Q;
              
end

I = I_sum/Niter;
Rmin = Rmin_sum/Niter;
Q = Q_sum/Niter;

hold on, plot(Rmin,Q*1e3,'g--')
hold on, plot(I,Q*1e3,'g')

xlabel('Rate [bits/channel use]'), ylabel('Harvested Energy [mJ]')
legend('1x1 (R_{lb},Q_{NL})','1x1 (R_{ca},Q_{NL})','2x1 (R_{lb},Q_{NL})','2x1 (R_{ca},Q_{NL})',...
       '4x1 (R_{lb},Q_{NL})','4x1 (R_{ca},Q_{NL})')
grid on

%% PS-TS on the same plot
clear; clc;

%-------------- CHANNEL --------------%
simulations = 1; 
if (simulations==1)
    Niter = 1e3;
    nT = 1; % number of transmit antennas
    kapa = nT; % shape parameter 
    theta = 1; % scale parameter 
    h = gamrnd(kapa,theta,1,Niter);
else
    Niter = 1;
    h = 1;
end

%-------------- CONSTANTS --------------%
% path loss
at = 0.5; % aperture of trasmit antenna [m]
ar = 0.01; % aperture of receive antenna [m]
fc = 2.45e9; % operating frequency
d = 12; % distance in [m] between BS and the user
L = 1-exp(-at*ar/d^2/(3e8/fc)^2); % path loss factor

% EH model: 1 for linear EH and 2 for nonlinear EH
model = 2;
if (model==1)
    z = 0.78; % enenrgy conversion efficiency
else
    A = 150; B = 0.014; Ps = 0.024;
end

% SNR at the receiver
P = 2*L; % [W]
snr_dB = 10; % [dB]
snr = 10^(snr_dB/10);
s = P/snr;
r = 3; % r = A/P

% average transmit power
if (r>2)
    E = P;
else
    E = r*P/2;
end

%-------------- RE REGION --------------%
sf = 0:1e-2:1; % splitting factor

IPS_sum = zeros(1,length(sf));
QPS_sum = zeros(1,length(sf));
ITS_sum = zeros(1,length(sf));
QTS_sum = zeros(1,length(sf));

for n = 1:length(h)    
    Nmp = max(floor(h(n)*r*snr/2.5),2);
    
    IPS = zeros(1,length(sf));
    QPS = zeros(1,length(sf));
    ITS = zeros(1,length(sf));
    QTS = zeros(1,length(sf));
    
    %-------------- PS --------------%
    for i = 1:length(sf)             
        sn = s*sqrt(1+1/(1-sf(i))^2); % standard deviation for PS 
        IPS(i) = Descrete_Mutual_Information(Nmp,h(n)*P,r,sn);
        if (model==1)
            QPS(i) = z*sf(i)*h(n)*E;
        else
            QPS(i)=(Ps/exp(A*B))*( (1+exp(A*B))/(1+exp(-A*(sf(i)*E*h(n)-B)))-1 );
        end

    end
    
    IPS_sum = IPS_sum + IPS;
    QPS_sum = QPS_sum + QPS;
    
    %-------------- TS --------------%
    for i = 1:length(sf)            
        sn = s*sqrt(2); % standard deviation
        ITS(i) = (1-sf(i))*Descrete_Mutual_Information(Nmp,h(n)*P,r,sn);
        if (model==1)
            QTS(i) = z*sf(i)*h(n)*E;
        else
            QTS(i)= sf(i)*(Ps/exp(A*B))*( (1+exp(A*B))/(1+exp(-A*(E*h(n)-B)))-1 );
        end
    end
    
    ITS_sum = ITS_sum + ITS;
    QTS_sum = QTS_sum + QTS;
              
end

IPS = IPS_sum/Niter;
QPS = QPS_sum/Niter;
IPS(end) = 0;

ITS = ITS_sum/Niter;
QTS = QTS_sum/Niter;

subplot(2,1,2)
plot(IPS,QPS*1e3,'b',ITS,QTS*1e3,'r')
xlabel('Rate [bits/channel use]'), ylabel('Harvested Energy [mJ]')
legend('PS scheme','TS scheme')
% xlim([0 0.7])
%max(max(IPS),max(ITS))

