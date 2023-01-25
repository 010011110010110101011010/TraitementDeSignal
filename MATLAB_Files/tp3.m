clear all 
close all 
clc 

%%
%Question 1:
load("ecg.mat");

fe = 500 ;
te = 1/fe;
N = length(ecg)
t = (0:N-1)*te;
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*(fe/N);


%% "Suppression du bruit provoqué par les mouvements du corps"

%Question 2
%  plot(t,ecg)
%  xlim([0.5 1.5])

%Question 3 : 
transF = fft(ecg) ;
% plot(f,transF)

%question 4 :
%generation d'un fltre pass_haut 
pass_haut = ones(size(ecg));
fc = 0.5;
index_fc = ceil((fc*N)/fe);
pass_haut(1:index_fc) = 0;
pass_haut(N-index_fc+1:N)=0;

% filtrage 
 
filtrage = pass_haut.*transF;
ecg1 = ifft(filtrage,"symmetric");

% apres le filtre pass_haut

%  plot(t,ecg1)
%  xlim([0.5 1.5])


%% "Suppression des interférences des lignes électriques 50Hz" 

%Question 6 : 
% generation d'u fitre Notch ideal

Notch_ideal = ones(size(ecg1));
fc_bande = 50 ; 
index_fc = ceil((fc_bande*N)/fe)+1;
Notch_ideal(index_fc) = 0 ;
Notch_ideal(N-index_fc+1) = 0 ;

tranFecg1 = fft(ecg1) ;
filter = Notch_ideal.*tranFecg1 ;
ecg2 = ifft(filter,"symmetric") ;
% hold on 
%  subplot(211)
%  plot(t,ecg2)
%   
%   subplot(212)
%  plot(t,ecg-ecg2)
%   plot(t,ecg1)
%  

% xlim([0.5 1.5])


%% "Amélioration du rapport signal sur bruit"

tranFecg2 = fft(ecg2) ;
% fshiftecg2 = fftshift(abs(2*tranFecg2/N));
% plot(fshift,fshiftecg2)

% Conception du filtre 

%Question 7 : 

pass_bas = zeros(size(ecg));
fc3 = 37 ;
index_fc3 = ceil((fc3*N)/fe);
pass_bas(1:index_fc3) = 1;
pass_bas(N-index_fc3+1:N)=1;

filter_bas = pass_bas.*tranFecg2 ;
ecg3 = ifft(filter_bas,"symmetric") ;

%Question 8 : 
subplot(211)
plot(t,ecg)
subplot(212)

plot(t,ecg-ecg3) 
%% /*"Identification de la fréquence cardiaque avec la fonction d’autocorrélation "*/

%Question 9 : 
[acf,lags] = xcorr(ecg3,ecg3);
stem(lags/fe,acf)
% Trouver le décalage correspondant à la fréquence cardiaque

%Question 10 : 
[maxval, maxidx] = max(acf);
local_maxima = find(acf(maxidx:end) < maxval, 1) + maxidx -1  ;

heart_rate = (fe / lags(local_maxima))/10;

