% Homework Exercise 3: 
% Validity check of Giannakis' formula
% Dafni Nikolaidou   10546

clc;
clear;
close all;

N = 2048;     % number of samples
q = 5;        % order of MA system

% MA process coefficients
b = [1.0, 0.93, 0.85, 0.72, 0.59, -0.1];

% Signal X construction
%[x,v] = X_Signal(N,b,q);  % uncomment to create a different X

load('X_signal.mat')

%% Q 1: Justify the non-Gaussian character of v (skewness calculation)
mu_v = mean(v);
std_v = std(v);
skewness_v = sum((v-mu_v).^3)/((N-1)*std_v^3);
fprintf('The skewness of v[k] is: %.4f\n', skewness_v);
% Not equal to zero thus v is non-Gaussian

%% Q 2:
% Estimation and plots of 3rd-order cumulants of x (indirect method K=32,
% M=64, L=20)
K = 32;
M = 64;
L = 20;

[~,~,cum3_x,~] = bisp3cum(x,M,L,'n','u');

% Plots of 3rd order cumulants
subplot(1,2,2)
axis=-L:L;
surf(axis,axis,cum3_x);
title('3rd order cumulants of x[k]');
xlabel('\tau_1'); ylabel('\tau_2');

subplot(1,2,1)
contour(axis,axis,cum3_x);
title('3rd order cumulants of x[k]');
xlabel('\tau_1'); ylabel('\tau_2');


%% Q 3:
% Estimate the impulse response of the MA system from Giannakis's formula
h = GiannakisFormula(cum3_x,q,L);
fprintf("Estimated impulse response of the MA system - Giannakis' formula" + ...
    " (order q=%d):\n",q);
disp(h);

%% Q 4:
% Estimation of h[k] from Giannakis's formula for sub-estimation and
% sup-estimation of the order q

% Sub-estimation of the order q
q_sub = q - 2;
h_sub = GiannakisFormula(cum3_x,q_sub,L);
fprintf("Estimated impulse response - Giannakis' formula" + ...
    "(q=%d):\n",q_sub);
disp(h_sub);

% Sup-estimation of the order q
q_sup = q + 3;
h_sup = GiannakisFormula(cum3_x,q_sup,L);
fprintf("Estimated impulse response - Giannakis' formula" + ...
    "(q=%d):\n",q_sup);
disp(h_sup);

%% Q 5:
% Estimate MA-q system output, plot original and estimated x[k]
% and find the NRMSE
x_est = conv(v,h,'same')';

figure();
plot(1:N,x);
hold on;
plot(1:N,x_est);
title('$x[k]$ and $\hat{x}[k]$ using Giannakis formula and $q=5$', 'Interpreter', 'Latex');
legend('Original signal','Estimated signal');

nrmse = NRMSE_calc(x,x_est,N);
fprintf('NRMSE for q = %d: %.4f\n',q,nrmse);

%% Q 6:
% Repeat Q 5 for h_sub and h_sup

% Sub-estimation of order q
x_est_sub = conv(v,h_sub,'same')';

figure();
plot(1:N,x);
hold on;
plot(1:N,x_est_sub);
title('$x[k]$ and $\hat{x}[k]$ using Giannakis formula and $q=3$', 'Interpreter', 'Latex');
legend('Original signal','Estimated signal');

nrmse_sub = NRMSE_calc(x,x_est_sub,N);
fprintf('NRMSE for q = %d: %.4f\n',q_sub,nrmse_sub);

% Sup-estimation of order q
x_est_sup = conv(v,h_sup,'same')';

figure();
plot(1:N,x);
hold on;
plot(1:N,x_est_sup);
title('x[k] and $\hat{x}[k]$ using Giannakis formula and $q=8$', 'Interpreter', 'Latex');
legend('Original signal','Estimated signal');

nrmse_sup = NRMSE_calc(x,x_est_sup,N);
fprintf('NRMSE for q = %d: %.4f\n',q_sup,nrmse_sup);

%% Q 7:
% Repeat Q 2,3,5 if we add AWGN at the output,
% producing a variation in SNR
snr = 30:-5:-5;
n = length(snr);
nrmse_awgn = zeros(1,n);

for i=1:n
    % output with noise
    y = awgn(x,snr(i),'measured');

    % 3rd order cumulants
    [~,~,cum3_y,~] = bisp3cum(y,M,L,'n','u');

    % impulse response - Giannakis' formula - q=5
    h_awgn = GiannakisFormula(cum3_y,q,L);

    % Output and NRMSE
    y_est = conv(v,h_awgn,'same')';
    nrmse_awgn(i) = NRMSE_calc(y,y_est,N);
end

% NRMSE - SNR plot
figure();
plot(snr,nrmse_awgn,'-o');
xlabel('SNR [dB]'); ylabel('NRMSE');
title('NRMSE in the output estimation vs SNR');

%% Q 8:
% Repeat the process 50 times and use mean values
r = 50;
NRMSE = zeros(1,r);
NRMSE_awgn = zeros(r,n);
NRMSE_sub = zeros(1,r); NRMSE_sup = zeros(1,r);

for i=1:r
    [x_i,v_i] = X_Signal(N,b,q);
    %3rd order cumulants
    [~,~,cum3_xi,~] = bisp3cum(x_i,M,L,'n','u');

    % Estimate impulse response from Giannakis' formula
    h_i = GiannakisFormula(cum3_xi,q,L);
    h_i_sub = GiannakisFormula(cum3_xi,q_sub,L);
    h_i_sup = GiannakisFormula(cum3_xi,q_sup,L);

    % Estimate output and NRMSE of each case
    x_est_i = conv(v_i,h_i,'same')';
    NRMSE(i) = NRMSE_calc(x_i,x_est_i,N);
    x_est_i = conv(v_i,h_i_sub,'same')';
    NRMSE_sub(i) = NRMSE_calc(x_i,x_est_i,N);
    x_est_i = conv(v_i,h_i_sup,'same')';
    NRMSE_sup(i) = NRMSE_calc(x_i,x_est_i,N);
    for j=1:n
        % output with noise
        y = awgn(x_i,snr(j),'measured');

        % 3rd order cumulants
        [~,~,cum3_y,~] = bisp3cum(y,M,L,'n','u');

        % impulse response - Giannakis' formula - q=5
        h_awgn = GiannakisFormula(cum3_y,q,L);

        % Estimate the output and find NRMSE
        y_est = conv(v_i,h_awgn,'same')';
        NRMSE_awgn(i,j) = NRMSE_calc(y,y_est,N);
    end
end
% Mean values of NRMSE
meanNRMSE = mean(NRMSE);
stdNRMSE = std(NRMSE);
meanNRMSE_sub = mean(NRMSE_sub);
stdNRMSE_sub = std(NRMSE_sub);
meanNRMSE_sup = mean(NRMSE_sup);
stdNRMSE_sup = std(NRMSE_sup);
meanNRMSE_awgn = mean(NRMSE_awgn);

fprintf(['\nMean values and standard deviation of NRMSE for %d ' ...
    'iterations:\n'],r);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q,meanNRMSE,stdNRMSE);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q_sub,meanNRMSE_sub, ...
    stdNRMSE_sub);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q_sup,meanNRMSE_sup, ...
    stdNRMSE_sup);

% 95% confidence interval - Bootstrap
alpha = 0.05;
B = 1000;     % bootstrap samples
boot_meanNRMSE_awgn = bootstrp(B, @mean, NRMSE_awgn);
% low and upper limits of confidence interval
low_lim = floor((B+1)*alpha/2);
up_lim = B+1-low_lim;

boot_meanNRMSE_awgn_sort = sort(boot_meanNRMSE_awgn);
bci_meanNRMSE_awgn = zeros(n,2);
for i = 1:n
    bci_meanNRMSE_awgn(i,1) = boot_meanNRMSE_awgn_sort(low_lim,i);
    bci_meanNRMSE_awgn(i,2) = boot_meanNRMSE_awgn_sort(up_lim,i);
end

% Plot mean NRMSE vs SNR
figure();
plot(snr,meanNRMSE_awgn,'-o');
hold on;
plot(snr,bci_meanNRMSE_awgn(:,1),'-.','Color','#A2142F');
plot(snr,bci_meanNRMSE_awgn(:,2),'-.','Color','#A2142F');
xlabel('SNR [dB]'); ylabel('NRMSE');
title('mean NRMSE in output estimation vs SNR');
legend('Mean value of NRMSE','95% confidence interval');