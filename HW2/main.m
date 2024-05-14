% Homework Exercise 2: 
% Bispectrum estimation using the indirect and the direct method
% Dafni Nikolaidou   10546

clc;
clear;
close all;

%% Q 1: Construction of X[k]
lambda = [0.12, 0.30, 0.42, 0.19, 0.17, 0.36];
omega = 2*pi*lambda;

% low and upper limits of the uniform distribution
low_lim = 0;
up_lim = 2*pi;
phi = (up_lim-low_lim)*rand(1,6)+low_lim;
phi(3) = phi(1) + phi(2);
phi(6) = phi(4) + phi(5);

N = 8192; % data length

% data
k = (0:N-1)';
X = sum(cos(k*omega+phi),2);
figure();
plot(X);
title('Real discrete process X[k]');

%% Q 2: Estimation of power spectrum and autocorrelation (128 max shiftings) 
muX = mean(X);
maxlag = 128;
corrX = xcorr(X,maxlag);
covariance = corrX-muX^2;

% Fourier transform
fs = 1;                        % sampling frequency
fftX = fft(covariance);
k = length(fftX);
f = (0:k-1)*(fs/k);            % frequency range
pow_spect = abs(fftX);
figure();
stem(f,pow_spect,'.');
xlabel('Frequency');
ylabel('Magnitude');
title('Power Spectrum');

%% Q 3a: 
% Estimation of bispectrum via the indirect method with K=32, M=256
% and L=64 max shiftings for the third order cumulants
K = 32;
M = 256;
L = 64;

% Bispectrum estimation - bispeci function (HOSA toolbox)
% i. Rectangular window
X = reshape(X,[M,K]);
figure();
bispecInd1 = bispeci(X,L,M,0,'unbiased',128,1);

% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% ii. Parzen window
figure();
bispecInd2 = bispeci(X,L,M,0,'unbiased',128,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

%% Q 3b: 
% Estimation of bispectrum via the direct method with K=32, M=256 and J=0
J = 0;
D = 2*J+1;
figure();
bispecDir1 = bispecd(X,M,D,M,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimated via the direct method');
legend('Bispectrum','Principal Region');

%% Q 7a: Repeat steps 3a and 3b for different segment lengths

%% i. 
% Estimation of bispectrum via the indirect method with K=16, M=512
% and L=64 max shiftings for the third order cumulants 
K = 16;
M = 512;

X = reshape(X,[M,K]);
% Rectangular window
figure();
bispecInd3 = bispeci(X,L,M,0,'unbiased',128,1);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% Parzen window
figure();
bispecInd4 = bispeci(X,L,M,0,'unbiased',128,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

% Estimation of bispectrum using the direct method with K=16, M=512 and J=0
figure();
bispecDir2 = bispecd(X,M,D,M,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the direct method');
legend('Bispectrum','Principal Region');

%% ii. 
% Estimation of bispectrum via the indirect method with K=64, M=128
% and L=64 max shiftings for the third order cumulants 
K = 64;
M = 128;

X = reshape(X,[M,K]);
% Rectangular window
figure();
bispecInd5 = bispeci(X,L,M,0,'unbiased',128,1);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% Parzen window
figure();
bispecInd6 = bispeci(X,L,M,0,'unbiased',128,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

% Estimation of bispectrum using the direct method with K=64, M=128 and J=0
figure();
bispecDir3 = bispecd(X,M,D,M,0);
% plotting the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Bispectrum estimation via the direct method');
legend('Bispectrum','Principal Region');

%% Q 7b:
% Repeat 3a and 3b for 50 realizations of the X[k] and compare
% the mean values of power spectrum and bispectrum
K = 32;
M = 256;
R = 50;
meanC2 = zeros(length(fftX),1);
meanC3Ind1 = zeros(M,M);
meanC3Ind2 = zeros(M,M);
meanC3Dir = zeros(M,M);

figure();
set(0,'DefaultFigureVisible','off');
for i=1:R
    % Construction of X process
    phi = (up_lim-low_lim)*rand(1,6)+low_lim;
    phi(3) = phi(1) + phi(2);
    phi(6) = phi(4) + phi(5);
    % data
    k = (1:N)';
    X = sum(cos(k*omega+phi),2);

    % Estimate autocorrelation (128 max shiftings) and power spectrum
    muX = mean(X);
    corrX = xcorr(X,maxlag);
    covariance = corrX-muX^2;
    % Fourier transform
    fftX = fft(covariance);
    pow_spect = abs(fftX);
    meanC2 = meanC2 + pow_spect;

    % Bispectrum estimation using indirect method
    X = reshape(X,[M,K]);
    % i. Rectangular window
    bispecInd = bispeci(X,L,M,0,'unbiased',128,1);
    meanC3Ind1 = meanC3Ind1 + bispecInd;
    % ii. Parzen window
    bispecInd = bispeci(X,L,M,0,'unbiased',128,0);
    meanC3Ind2 = meanC3Ind2 + bispecInd;

    % Estimation of bispectrum using direct method
    bispecDir = bispecd(X,M,D,M,0);
    meanC3Dir = meanC3Dir + bispecDir;
end

meanC2 = meanC2/R;
meanC3Ind1 = meanC3Ind1/R;
meanC3Ind2 = meanC3Ind2/R;
meanC3Dir = meanC3Dir/R;

nfft = 256;
if rem(nfft,2) == 0
    waxis = (-nfft/2:(nfft/2-1))/nfft;
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)/nfft;
end

% Plot the power spectrum
fs = 1;
k = length(meanC2);
f = (0:k-1)*(fs/k);
set(0,'DefaultFigureVisible','on');
set(gcf,'Name','');
plot(f,meanC2);
xlabel('Frequency');
ylabel('Magnitude');
title('Mean value of power spectrum estimations');

% Plot the bispectrum (indirect method, rectangular window)
figure();
contour(waxis,waxis,abs(meanC3Ind1),4); grid on;
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title(['Mean estimation of bispectrum estimated via the indirect ' ...
    'method - Rectangular window']);
xlabel('f1');
ylabel('f2');

% Plot the bispectrum (indirect method, Parzen window)
figure();
contour(waxis,waxis,abs(meanC3Ind2),4); grid on;
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title(['Mean estimation of bispectrum estimated via the indirect ' ...
    'method - Parzen window']);
xlabel('f1');
ylabel('f2');

% Plot the bispectrum (direct method)
figure();
contour(waxis,waxis,abs(meanC3Dir),4); grid on;
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','r');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','r');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','r');                % f2=0
title('Mean estimation of bispectrum via the direct method');
xlabel('f1');
ylabel('f2');
