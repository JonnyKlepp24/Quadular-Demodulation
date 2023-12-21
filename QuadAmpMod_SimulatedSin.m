%% Quadrature Amplitude Modulation from Simulated Signal w/noise
% a = b * c
% c = b * ^-1[a] -> c = b/a

clear all; clc; close all;

%% Generation of multifrequency wave

N = 2048;           % Number of point
k = 0:N-1;          % Vector of number of points 
f_samp = 200000;    % Sample frequency 
tk = k/f_samp;      % X vector of 1/f_samp space

% Frequency of waves
w01 = 2*pi*10000;
w02 = 2*pi*30000;
w03 = 2*pi*50000;

phase = pi/6; % Phase
A = 1;        % Amplitude
c = 0.1;      % Offset

% Assmeble the signal
signalClean = A*sin(w01*tk) + A*sin(w02*tk) + A*sin(w03*tk);

% Add noise to signal
figure
signalNoise = awgn(signalClean,10,'measured');

% Plot clean signal vs noisy signal
tiledlayout(2,1)
nexttile
plot(tk',[signalClean' signalNoise'])
legend('Original Signal','Signal with AWGN')
nexttile
plot(tk',[signalClean' signalNoise'])
legend('Original Signal','Signal with AWGN')
xlim([1/w03 200/w03])

%% Signal extraction

% matrix containing the desired frequency component of the system
E1 = [sin(w01*tk)', cos(w01*tk)', ones(1, length(tk))']; 
E2 = [sin(w02*tk)', cos(w02*tk)', ones(1, length(tk))'];
E3 = [sin(w03*tk)', cos(w03*tk)', ones(1, length(tk))'];
Etot = [E1 E2 E3];

% Extract all parts of phi
phi_tot = pinv(Etot)*signalNoise';
% inverse of E * s^t (pinv)
%phi_tot =  E\signal';

% recreate each signal at its freqeuncy
demodF1 = Etot(:,1:3) * phi_tot(1:3,:);   % signal @ frequency 1
demodF2 = Etot(:,4:6) * phi_tot(4:6,:);   % signal @ frequency 2
demodF3 = Etot(:,7:9) * phi_tot(7:9,:);   % signal @ frequency 3

% combine all the different signal
signalSummedDemod = demodF1 + demodF2 + demodF3;

%% Seperate phi characteristics
alpha1 = phi_tot(1);
beta1 = phi_tot(2);
Cout1 = phi_tot(3);

alpha2 = phi_tot(4);
beta2 = phi_tot(5);
Cout2 = phi_tot(6);

alpha3 = phi_tot(7);
beta3 = phi_tot(8);
Cout3 = phi_tot(9);

%% Evaluation
% Error
dif = signalClean' - signalSummedDemod;
err = ((dif)')*(dif);

% frequency 1
offset1 = mean(demodF1)
ampVp1 = sqrt(alpha1^2 + beta1^2)
phase1 = asind(beta1/ampVp1)

% frequency 2
offset2 = mean(demodF2)
ampVp2 = sqrt(alpha2^2 + beta2^2)
phase2 = asind(beta2/ampVp2)

% frequency 3
offset3 = mean(demodF3)
ampVp3 = sqrt(alpha3^2 + beta3^2)
phase3 = asind(beta3/ampVp3)

%% PLOTS

figure
subplot(2,1,1);
plot(tk, signalNoise, 'color', 'b')
title("Measured Signal (noisy)")
xlabel("Time (s)");
ylabel("Voltage (V)");
subplot(2,1,2);
plot(tk, signalNoise,'.-','color', 'b')
xlim([1/w03 200/w03])
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
subplot(2,1,1);
plot(tk, signalSummedDemod, 'color', 'b')
title("Summed Demodulated Signal")
xlabel("Time (s)");
ylabel("Voltage (V)");
subplot(2,1,2);
plot(tk, signalSummedDemod,'.-','color', 'b')
xlim([1/w03 200/w03])
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
plot(tk, dif)
title("Measured - Summed Demodulated")
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
plot(tk,signalSummedDemod)
title("Summed Multifrequncy Signal")
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
subplot(2,1,1);
hold on;
plot(tk, demodF1, 'color', 'b')
title("E1 Modulated Signal")
xlabel("Time (s)");
ylabel("Voltage (V)");
subplot(2,1,2);
hold on;
plot(tk, demodF1,'.-','color', 'b')
xlim([1/w03 50/w03])
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
subplot(2,1,1);
hold on;
plot(tk, demodF2, 'color', 'b')
title("E2 Modulated Signal")
xlabel("Time (s)");
ylabel("Voltage (V)");
subplot(2,1,2);
hold on;
plot(tk, demodF2,'.-','color', 'b')
xlim([1/w03 50/w03])
xlabel("Time (s)");
ylabel("Voltage (V)");

figure
subplot(2,1,1);
hold on;
plot(tk, demodF3, 'color', 'b')
title("E3 Modulated Signal")
xlabel("Time (s)");
ylabel("Voltage (V)");
subplot(2,1,2);
hold on;
plot(tk, demodF3,'.-','color', 'b')
xlim([1/w03 50/w03])
xlabel("Time (s)");
ylabel("Voltage (V)");
