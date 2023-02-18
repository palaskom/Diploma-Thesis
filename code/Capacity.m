%% plot v(r(mu))
clear;clc
r = 1:0.01:8;
v = zeros(1,length(r));
for i = 1:length(r)
    [~,v(i)] = CapacityLowerBound(r(i),1,1);
end
plot(r,v), xlabel('$r = A/P$','Interpreter','latex'), ylabel('$v$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)

%% Plot Capacity Lower and Upper Bounds (in [nats/channel use])

%----------------------------------------------%
%        change all log2() -> log()            %
%----------------------------------------------%

clear;clc

sigma = 1;
P = 10^(0):0.05:10^(2);

% upper and lower bounds
Rmin = zeros(1,length(P));
Rmax = zeros(1,length(P));

r = 2; % r = A/P
for i = 1:length(P)
    Rmin(i) = CapacityLowerBound(r,P(i),sigma);
    Rmax(i) = CapacityUpperBound(r,P(i),sigma);
end
subplot(2,1,1)
plot(10*log10(P),Rmax,'k',10*log10(P),Rmin,'r')
legend('Upper bound','Lower bound'), grid on
xlabel('SNR [dB]'), ylabel('Mutual Information [nats/channel use]')

r = 6; % r = A/P
for i = 1:length(P)
    Rmin(i) = CapacityLowerBound(r,P(i),sigma);
    Rmax(i) = CapacityUpperBound(r,P(i),sigma);
end
subplot(2,1,2)
plot(10*log10(P),Rmax,'k',10*log10(P),Rmin,'r')
legend('Upper bound','Lower bound'), grid on
xlabel('SNR [dB]'), ylabel('Mutual Information [nats/channel use]')


%% Plot Mutual Information along with Lower and Upper Bounds (in [nats/channel use])

%----------------------------------------------%
%        change all log2() -> log()            %
%----------------------------------------------%

clear;clc

sigma = 1;

P = 10^(-0.6):0.05:10^(1.4);
r = 3; % r = A/P
Nmp = 2:2:8; % number of mass points

% mutual information of the discrete distribution
for k = 1:length(Nmp)
    I = zeros(1,length(P));
    for i = 1:length(P)
        I(i) = Descrete_Mutual_Information(Nmp(k),P(i),r,sigma);    
    end   
    hold on, plot(10*log10(P),I)
end

% upper and lower bounds
Rmin = zeros(1,length(P));
Rmax = zeros(1,length(P));
for i = 1:length(P)
    Rmin(i) = CapacityLowerBound(r,P(i),sigma);
    Rmax(i) = CapacityUpperBound(r,P(i),sigma);
end
hold on, plot(10*log10(P),Rmax,'k--')    
hold on, plot(10*log10(P),Rmin,'r--')
% xlim([0,15]) 
% xticks(0:2.5:15)
legend('R_{ca} (2 mass points)','R_{ca} (4 mass points)','R_{ca} (6 mass points)','R_{ca} (8 mass points)',...
       'R_{ub} (upper bound)','R_{lb} (lower bound)')
xlabel('SNR [dB]'), ylabel('Mutual Information [nats/channel use]')
grid on

%% Plot Mutual Information vs SNR for different values of r (in [bits/channel use])

clear;clc

sigma = 1;
snr_dB = 5:0.1/2:15; %[dB]
snr = 10.^(snr_dB/10);
P = snr/sigma;
r = 10; % r = A/P

I = zeros(1,length(snr));
for i = 1:length(I)
    if (snr_dB(i)<=4) 
        c = 2.2;
    else
        c = 2.5;
    end
    Nmp = max(2,floor(r*snr(i)/c)); 
    I(i) = Descrete_Mutual_Information(Nmp,P(i),r,sigma); 
end
hold on, plot(snr_dB,I)
xlabel('SNR [dB]'), ylabel('Mutual Information [bits/channel use]')
grid on
legend({'$r = 4$','$r = 6$','$r = 10$'},'Interpreter','latex')