clc
clear all 
close all
                                        %   Filtrage et diagramme de Bode

te=0.0005;
fe=1/te;

f1=500;
f2=400;
f3=50;

t=0:te:5;
N=length(t);

f =  (0:N-1)*(fe/N);
freqshift=linspace(-fe/2,fe/2,N);

X=sin(2*pi*f1*t)  + sin(2*pi*f2*t)  + sin(2*pi*f3*t);
x_Spectrum=abs(fft(X))/N;

figure
    subplot(2,1,1)
        plot(t,X)
            xlabel('Temps (s)');ylabel('Amplitude')
            title('Signal')
    subplot(2,1,2)
        plot(freqshift,fftshift(x_Spectrum)) 
            xlabel('Amplitude ');ylabel('Frequence (Hz)')
            title('Spectre du signal')  
            
% Compatiblement avec la formule de signal X son spectre illustre que ce dernier est constitué de trois composantes fréquentielles d'amplitudes égales       
%%  Diagramme de Bode  de la transmittance complexe  H(f)

k=1;
Wc=50;
w = 2*pi*f;

H=(k*j*w/Wc)./(1+j*w/Wc);

% figure
% H = tf([Wc],[1 Wc]);
% bode(H)

amplitude = 20*log(abs(H));
phase = angle(H);
figure
    semilogx(amplitude)
    xlabel('Frequence (rad/s)');ylabel('Amplitude')
    title('Amplitude de la transmittance complexe H(f)')
        
            
%%  Representation  des amplitudes des trois transmittances complexes
wc1 = 2*pi*100;
wc2 = 2*pi*250;
wc3 = 2*pi*500;

H1 = (k*j*w/wc1) ./ (1 + 1j*w/wc1);
H2 = (k*j*w/wc2) ./ (1 + 1j*w/wc2);
H3 = (k*j*w/wc3) ./ (1 + 1j*w/wc3);

Amp1=20*log(abs(H1));
Amp2=20*log(abs(H2));
Amp3=20*log(abs(H3));

Ph1=angle(H1);
Ph2=angle(H2);
Ph3=angle(H3);

figure
    subplot(2,1,1)
        semilogx(f,Amp1,t,Amp2,t,Amp3,"Linewidth",1.5)
            xlabel('Amplitude ');ylabel('Frequence (Hz)')
            title('Module de la fonction H(t) ')
            legend(" fc=100","fc=250","fc=500")
            grid on
    subplot(2,1,2)
        semilogx(f,Ph1,f,Ph2,f,Ph3)
            xlabel('Phase ');ylabel('Frequence (Hz)')
            title('Modules des fonctions H(t) ')
            legend(" fc=100","fc=250","fc=500")
            grid on
            
% Appliquer le filtre au signal x(t)

% filtre avec fc=100
X_Filtered1 = ifft(H1.*x_Spectrum,'symmetric');

% filtre avec fc=250
X_Filtered2 = ifft(H2.*x_Spectrum,'symmetric');

% filtre avec fc=500
X_Filtered3 = ifft(H3.*x_Spectrum,'symmetric');

figure
    subplot(2,2,1)
        plot(freqshift, fftshift(x_Spectrum))
            xlabel('Frequence (Hz) ');ylabel('Amplitude')
            title('Spectre du signal initial ')    
    subplot(2,2,2)
        plot(freqshift, fftshift(abs(fft(X_Filtered1))))
            xlabel('Frequence (Hz) ');ylabel('Amplitude')
            title('Spectre du signal filtré  fc=100')
    subplot(2,2,3)
         plot(freqshift, fftshift(abs(fft(X_Filtered2))))
            xlabel('Frequence (Hz) ');ylabel('Amplitude')
            title('Spectre du signal filtré  fc=250')
    subplot(2,2,4)
        plot(freqshift, fftshift(abs(fft(X_Filtered3))))
            xlabel('Frequence (Hz)');ylabel('Amplitude')
            title('Spectre du signal filtré  fc=500 ')

% Le filtre le plus efficace est celui ayant une fréquence de coupure de 250 hz 
%or au-delà de cette fréquence les autres composantes fréquentielles commencent à être infectées


%% 
figure
    subplot(2,1,1)
        plot(freqshift, X)
            xlabel('Frequence (Hz) ');ylabel('Amplitude')
            title('Spectre du signal initial ')    
    subplot(2,1,2)
        plot(freqshift, X_Filtered1)
            xlabel('Frequence (Hz) ');ylabel('Amplitude')
            title('Spectre du signal filtré  fc=100')


%%   Dé-bruitage d'un signal sonore
clc
%fc= f tel que g    ain est de  1/sqr(2)
[Song,fs]=audioread('test.wav');
% sound(Song,fs)

% Puisque on est dit que le bruit est d'une haute frequence donc on va le
% filtrer en faisant recours a un filtre pass bas 


% Définissons la fréquence de coupure
%fc = 2000;

% Créons la fonction de transfert du filtre
%[B, A] = butter(4, fc/(fs/2), 'low');

% Appliquons le filtre au signal enregistré
%Filtered_Song = filter(B, A, Song);

%sound(Filtered_Song,fs)


% temps=0:1/fs:5;
% n=length(temps);
% 
% frequence =  (0:n-1)*(fs/n);
% W=2*pi*frequence;
% 
% W_cutoff =2000;
% 
% H_Passe_Bas=1./(1+j*W/W_cutoff);
% Song_Spectrum=fft(Song);
% 
% Filtered_Song=H_Passe_Bas.*fft(Song)
% sound(Filtered_Song,fs)