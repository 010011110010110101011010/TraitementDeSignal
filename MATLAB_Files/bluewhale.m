clear all
close all
clc

%%

%Q1 : chagement du fichier 

[x,fs]= audioread("bluewhale.au"); 
chant = x(2.45e4:3.10e4); 

%% 
% Definition des variables 

N = length(chant); 
ts = 1/fs; 
t = (0:N-1)*(10*ts); 
f = (0:floor(N/2))*(fs/N)/10;

%%

% Q2 : visualisation du l'ecout du signal 
sound(chant)
plot(t,chant);

%%

%Q3 : densite spectrale de signal

dsp_chant =  (abs(fft(chant)).^2)/N;

%  Q4 : Déterminer à partir du tracé, la fréquence fondamentale du gémissement de rorqual 

plot(f,dsp_chant(1:floor(N/2)+1))

% nous pouvons remarquer  que la ff est 16,9205 
