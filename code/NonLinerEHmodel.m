%% proposed nonlinear EH model

clear; clc

p1 = 0.005; p2 = 0.04; dp = p2-p1;
P = p1:dp/100:p2;

% nonlinear model
Ps = 0.024;
A = 150; B = 0.014; % @ 2.45 GHz
y = 1./(1+exp(-A*(P-B)));
c = 1/(1+exp(A*B));
Qnl = Ps*(y-c)/(1-c);

subplot(2,1,1) % plot nonlinear model
plot(P*10^3,Qnl*10^3)
xlabel('Input RF power [mW]'), ylabel('Harvested Power [mW]')

% linear fitting to the above nonlinear model for Pin [5,25]mW
p1 = 0.005; p2 = 0.025; dp = p2-p1;
P = p1:dp/100:p2;

Ps = 0.024;
A = 150; B = 0.014; % @ 2.45 GHz
y = 1./(1+exp(-A*(P-B)));
c = 1/(1+exp(A*B));
Qnl = Ps*(y-c)/(1-c);

F = @(a,P) a*P;
a = lsqcurvefit(F,0.8,P,Qnl);

% linear model
z = a;
Ql = z*P;

subplot(2,1,2)
plot(P*10^3,Qnl*10^3,'b',P*10^3,Ql*10^3,'r')
% figure, plot(P*10^3,Qnl*10^3), hold on, plot(P*10^3,Ql*10^3,'r')
xlabel('Input RF power [mW]'), ylabel('Harvested Power [mW]')
legend('nonlinear EH model','linear EH model')

% z = Q./P;
% figure, plot(P*10^3,z*100)
% xlabel('Input RF power (mW)'), ylabel('Efficiency (%)')
% xlim([P(1)*10^3 P(end)*10^3])