%% PS/TS scheme - plot Po(d,nT)

%-----------------------------------------------------------%
% Plot of Joint Harvested Energy-Rate Outage Probability as  
% a function of the distance between the BS and the user and
% the number of transmit antennas
%-----------------------------------------------------------%

clear;clc

P = 2; % average power constraint in [W]
r = 3; % peak-to-average power constraint A/P
rth = 0.1; % rate threshold in [bits/channel use]
qth = 0.1e-3; % energy threshold in [J]

% noise power in [W]: 
% sigma = srec = sadc = L*P/snr, where L is computed @d=12m and snr=10dB
srec = 4.4395e-04;  
sadc = srec; 

d = 5:0.05:15; % distance between BS and the user in [m]
Po_PS = zeros(1,length(d));
Po_TS = zeros(1,length(d));
fPS = zeros(1,length(d));
fTS = zeros(1,length(d));

nT = [1 2 4]; % number of trasmit antennas
C = {'r','b','g'}; % colors used for plot

for j = 1:length(nT)
    a = nT(j); 
    b = 1;
    for i = 1:length(d)
        [Po_PS(i),fPS(i)] = OutProbPs(P, r, d(i), rth, qth, sadc, srec, a, b, 2);
        [Po_TS(i),fTS(i)] = OutProbTs(P, r, d(i), rth, qth, sadc, srec, a, b, 2);
    end
    
    semilogy(d,Po_PS,C{j}), hold on
    semilogy(d,Po_TS,'LineStyle','--','Color',C{j})
    
%     % splitting factor remains the same when the distance between BS and the user changes
%     hold on, plot(d,fPS,C{j})
%     hold on, plot(d,fTS,'LineStyle','--','Color',C{j})
end

xlabel('Distance between BS and the user [m]')
ylabel('Joint Outage Probability')
legend('1x1 (PS scheme)','1x1 (TS scheme)','2x1 (PS scheme)','2x1 (TS scheme)','4x1 (PS scheme)','4x1 (TS scheme)')
grid on


%% PS/TS scheme - plot Po(P,nT)

%-----------------------------------------------------------%
% Plot of Joint Harvested Energy-Rate Outage Probability as 
% a function of the average power P and the number of transmit 
% antennas.
% Note: P is the average power constraint and not the average
% transmit power. The average transmit power is E and is used
% in OutProbPS/TS.m functions. If r>2, E=P. If r<=2, E=r*P/2
%-----------------------------------------------------------%

clear;clc

r = 3; % peak-to-average power constraint A/P
d = 12; % distance between BS and the user in [m]
rth = 0.1; % rate threshold in [bits/channel use]
qth = 0.1e-3; % energy threshold in [J]

% noise power in [W]: 
% sigma = srec = sadc = L*P/snr, where L is computed @d=12m and snr=10dB
srec = 4.4395e-04;  sadc = srec; 

P = 1:0.01:6; % average power constraint in [W]

fPS = zeros(1,length(P));
fTS = zeros(1,length(P));
Po_PS = zeros(1,length(P));
Po_TS = zeros(1,length(P));

nT = [1 2 4]; % number of trasmit antennas
C = {'r','b','g'};

for j = 1:length(nT)   
    a = nT(j); 
    b = 1;
    for i = 1:length(P)
        [Po_PS(i),fPS(i)] = OutProbPs(P(i), r, d, rth, qth, sadc, srec, a, b, 2);
        [Po_TS(i),fTS(i)] = OutProbTs(P(i), r, d, rth, qth, sadc, srec, a, b, 2);
    end   
    semilogy(P,Po_PS,C{j}), hold on
    semilogy(P,Po_TS,'LineStyle','--','Color',C{j})
%     % splitting factor remains the same when the average transmit power changes
%     hold on, plot(P,fPS,C{j})
%     hold on, plot(P,fTS,'LineStyle','--','Color',C{j})
    
end

% xticks([0.1 0.5:0.5:5])
grid on, xlim([P(1) P(end)])
xlabel('Average Transmit Power [W]'), ylabel('Joint Outage Probability')
legend('1x1 (PS scheme)','1x1 (TS scheme)','2x1 (PS scheme)','2x1 (TS scheme)','4x1 (PS scheme)','4x1 (TS scheme)')

%% PS/TS scheme - plot Po(LP/sigma,nT)

%-----------------------------------------------------------%
% Plot of Joint Harvested Energy-Rate Outage Probability as 
% a function of the average power P and the number of transmit 
% antennas.
% Note: P is the average power constraint and not the average
% transmit power. The average transmit power is E and is used
% in OutProbPS/TS.m functions. If r>2, E=P. If r<=2, E=r*P/2
%-----------------------------------------------------------%

clear;clc

r = 3; % peak-to-average power constraint A/P
d = 12; % distance between BS and the user in [m]
rth = 0.1; % rate threshold in [bits/channel use]
qth = 0.1e-3; % energy threshold in [J]

at = 0.5; % aperture of trasmit antenna [m]
ar = 0.01; % aperture of receive antenna [m]
fc = 2.45e9; % operating frequency
L = 1-exp(-at*ar/d^2/(3e8/fc)^2); % path loss factor

snr = 5:0.05:15;
s = 4.4395e-04;
P = (s/L)*10.^(snr/10);

fPS = zeros(1,length(P));
fTS = zeros(1,length(P));
Po_PS = zeros(1,length(P));
Po_TS = zeros(1,length(P));

nT = [1 2 4]; % number of trasmit antennas
C = {'r','b','g'};

for j = 1:length(nT)   
    a = nT(j); 
    b = 1;
    for i = 1:length(P)
        [Po_PS(i),fPS(i)] = OutProbPs(P(i), r, d, rth, qth, s, s, a, b, 2);
        [Po_TS(i),fTS(i)] = OutProbTs(P(i), r, d, rth, qth, s, s, a, b, 2);
    end   
    semilogy(snr,Po_PS,C{j}), hold on
    semilogy(snr,Po_TS,'LineStyle','--','Color',C{j})
%     % splitting factor remains the same when the average transmit power changes
%     hold on, plot(P,fPS,C{j})
%     hold on, plot(P,fTS,'LineStyle','--','Color',C{j})
    
end

% xticks([0.1 0.5:0.5:5])
% grid on, xlim([P(1) P(end)])
% yticks([10^(-1) 10^(-2) 10^(-2) 10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)])
% yticks([10 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9]*1e-9)
yticks([10 10^3 10^5 10^7 10^9]*1e-9)
ylabel('Joint Outage Probability')
legend('1x1 (PS scheme)','1x1 (TS scheme)','2x1 (PS scheme)','2x1 (TS scheme)','4x1 (PS scheme)','4x1 (TS scheme)')
% xlabel('$$\frac{LP}{\sigma} [dB]$$','fontsize',10,'Interpreter','Latex')
xlabel('LP/\sigma [dB]')
grid on

%% PS/TS scheme - plot Po(r,nT)
clear;clc

d = 12; % distance between BS and the user in [m]
P = 2; % average power constraint in [W]
rth = 0.1; % rate threshold in [bits/channel use]
qth = 0.1e-3; % energy threshold in [J]

% noise power in [W]: 
% sigma = srec = sadc = L*P/snr, where L is computed @d=12m and snr=10dB
srec = 4.4395e-04;  sadc = srec; 

r = 1:0.01:6; % peak-to-average power constraint A/P

fPS = zeros(1,length(r));
fTS = zeros(1,length(r));
Po_PS = zeros(1,length(r));
Po_TS = zeros(1,length(r));

C = {'r','b','g'};
nT = [1 2 4]; % number of trasmit antennas

for j = 1:length(nT)
    a = nT(j);
    b = 1;
    for i = 1:length(r)
        [Po_PS(i),fPS(i)] = OutProbPs(P, r(i), d, rth, qth, sadc, srec, a, b, 2);
        [Po_TS(i),fTS(i)] = OutProbTs(P, r(i), d, rth, qth, sadc, srec, a, b, 2);
    end    

%     semilogy(r,Po_PS,C{j}), hold on
%     semilogy(r,Po_TS,'LineStyle','--','Color',C{j})
    hold on, plot(r,fPS,'b')
    hold on, plot(r,fTS,'r')
end

% xlabel('Peak-to-Average Power Constraint r = A/P')
% ylabel('Joint Outage Probability')
% legend('1x1 (PS scheme)','1x1 (TS scheme)','2x1 (PS scheme)','2x1 (TS scheme)','4x1 (PS scheme)','4x1 (TS scheme)')
% grid on

xlabel('Peak-to-Average Power Ratio r = A/P')
ylabel('Optimal splitting factor')
legend('PS scheme','TS scheme')
grid on


%% Outage Probability - (rth,qth) 
clear;clc

d = 12; % distance between BS and the user in [m]
P = 2; % average power constraint in [W]
r = 3; % peak-to-average power constraint A/P

% noise power in [W]: 
% sigma = srec = sadc = L*P/snr, where L is computed @d=12m and snr=10dB
srec = 4.4395e-04;  sadc = srec; 

nT = 4; % number of trasmit antennas
a = nT; b = 1;

rth = linspace(0.05,0.6,200); % rate threshold in [bits/channel use] 0.01,1,100
qth = linspace(0.02,0.5,200)*1e-3; % energy threshold in [J] 0.05,0.5,100
[Rth,Qth] = meshgrid(rth,qth);
Po_PS = zeros(length(qth),length(rth)); % outage probability for PS scheme
Po_TS = zeros(length(qth),length(rth)); % outage probability for PS scheme

for i = 1:length(qth) % fixed row => fixed qth
    for j = 1:length(rth) % fixed row => fixed rth
        Po_PS(i,j) = OutProbPs(P, r, d, rth(j), qth(i), sadc, srec, a, b, 2);
        Po_TS(i,j) = OutProbTs(P, r, d, rth(j), qth(i), sadc, srec, a, b, 2);
    end
end

% Outage Probability for PS scheme - (rth,qth)
figure
pcolor(Rth,Qth*1e3,Po_PS), shading flat; colorbar;
xlabel('Rate threshold [bits/channel use]'), ylabel('Energy threshold [mJ]')
set(gca,'ColorScale','log')
xticks(0.05:0.05:0.6)
yticks(0.02:0.06:0.5)
% c = colorbar;
% c.Ticks = linspace(min(min(Po_PS)),max(max(Po_PS)),10);
% xlim([min(min(Rth)) max(max(Rth))])
% c.Ruler.Scale = 'log';
% c.Ruler.MinorTick = 'on';
% set(gca,'ZScale','log')

% Outage Probability for TS scheme - (rth,qth)
figure
pcolor(Rth,Qth*1e3,Po_TS), shading flat; colorbar;
xlabel('Rate threshold [bits/channel use]'), ylabel('Energy threshold [mJ]')
set(gca,'ColorScale','log')
xticks(0.05:0.05:0.6)
yticks(0.02:0.06:0.5)
% c = colorbar;
% c.Ticks = linspace(min(min(Po_TS)),max(max(Po_TS)),10);


%% Optimal PS/TS factor - (rth,qth) 
clear;clc

d = 12; % distance between BS and the user in [m]
P = 2; % average power constraint in [W]
r = 3; % peak-to-average power constraint A/P

% noise power in [W]: 
% sigma = srec = sadc = L*P/snr, where L is computed @d=12m and snr=10dB
srec = 4.4395e-04;  sadc = srec; 

nT = 1; % number of trasmit antennas
a = nT; b = 1;

rth = linspace(0.05,0.5,200); % rate threshold in [bits/channel use] 0.01,1,100
qth = linspace(0.05,0.5,200)*1e-3; % energy threshold in [J] 0.05,0.5,100
[Rth,Qth] = meshgrid(rth,qth);
fPS = zeros(length(qth),length(rth)); % PS splitting factor value on (rth,qth)
fTS = zeros(length(qth),length(rth)); % TS splitting factor value on (rth,qth)

for i = 1:length(qth) % fixed row => fixed qth
    for j = 1:length(rth) % fixed row => fixed rth
        [~,fPS(i,j)] = OutProbPs(P, r, d, rth(j), qth(i), sadc, srec, a, b, 2);
        [~,fTS(i,j)] = OutProbTs(P, r, d, rth(j), qth(i), sadc, srec, a, b, 2);
    end
end

% Optimal PS factor - (rth,qth)
figure
pcolor(Rth,Qth*1e3,fPS), shading flat; colorbar;
xlabel('Rate threshold [bits/channel use]'), ylabel('Energy threshold [mJ]')
c = colorbar;
c.Ticks = linspace(min(min(fPS)),max(max(fPS)),10);
xticks(0.05:0.05:0.5)

% Optimal TS factor - (rth,qth)
figure
pcolor(Rth,Qth*1e3,fTS), shading flat; colorbar;
xlabel('Rate threshold [bits/channel use]'), ylabel('Energy threshold [mJ]')
c = colorbar;
c.Ticks = linspace(min(min(fTS)),max(max(fTS)),10);
xticks(0.05:0.05:0.5)
