clear all
close all
clc


%%
%Definition des variables 

fe = 1e4; % la frequence d’échantillonnage
te = 1/fe; % le temps  d’échantillonnage
N = 5000; % le nombre d'echantillon 
t = (0:N-1)*te; %temps d'échantillonnage pour chaque échantillon
f = (0:N-1)*(fe/N);% fréquences pour chaque échantillon
fshift = (-N/2:N/2-1)*(fe/N); %un échantillonnage fréquentiel d'un signal qui a été décalé
%%

%definition du signal  x(t)

x = 1.2*cos(2*pi*440*t+1.2)+ 3*cos(2*pi*550*t)+ cos(2*pi*2500*t);

%%

%Q1 : tracage du signal dans le domaine temporelle : 

plot(t,x)


%%

%Q2 : Calculons le TFD du signal puis tracons son spectre d'amplitude : 

%la commande fft nous a permet de calculer le transformer de Fourier 
%qui est  le sectre 

spectre_x = fft(x);

% on sait que l'amplitude est le modele de TF pour l'obtenir on utlise
% la commande abs 
spectre_amplitude = abs(spectre_x) ;

% Nous pouvons remarquer que le sepectre est centarlise par N 
plot(f,spectre_amplitude) ;

%%

%Q3 : Utilisation de fftshift Pour mieux visualiser le contenu fréquentiel du signal

% avec la commande fftshift on a effectue un decalage circulaire avec lequel 
%on a obtenu un spectre censtarlise par 0 

plot(fshift , 2*fftshift(spectre_amplitude/N) )


%%

%Q4 : generation d'un signal xnoise 
% Ce bruit suit une loi normale de moyenne 0 et d’écart type 1
 noise = randn(size(x));
 plot(noise)
% ajoutons le bruit à notre signal 
 xnoise = x + noise ;

 
%%

%Q5 : écout du  signal et puis le signal bruité. 
% sound(noise)
%sound(xnoise)

%%

%Q6 : le spectre  du signal bruité
ynoise = fft(xnoise);
plot(fshift,fftshift(2*abs(ynoise)/N))

%%

% Conception du filtre 

pass_bas = zeros(size(x));
fc = 2000;
index_fc = ceil((fc*N)/fe);
 pass_bas(1:index_fc) = 1;
pass_bas(N-index_fc+1:N)=1;
% 
% 
plot(f,pass_bas,"linewidth",1.5)
% 
% 
% 
% % Filtrage
% 
 x_filter_freq = pass_bas.*spectre_x;
 x_filter_tmp = ifft(x_filter_freq,"symmetric");
% 
% 
  plot(fshift, fftshift(abs(fft(x_filter_tmp))));
% 




