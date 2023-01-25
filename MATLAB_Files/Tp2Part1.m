clear all
close all
clc
%% Jeux de mots
%Question 1:

[x,fe]=audioread('sound.au');

%Question 2:

%sound(x,fe)
%plot(x)

%Question 3:

%sound(x,fe*2)
%sound(x,fe/2)
%sound(x,fe)



%Question 4:

% stem(x)
%Après le traçage de la fonction x on peut arriver a voir que le morceau 
%Rien ne sert de est contenue dans la partie de l'index 2695 j'usqu'à 53180

%Question 5:


riennesertde=x(2695:(53180));
% stem(riennesertde)
% sound(riennesertde,fe)

%Question 6:
% 
 courir=x((53311):(69800));
%  sound(courir,fe)
% 
ilfaut=x((72130):(81243));
% sound(ilfaut,fe)
% 
partirapoint=x((86700):(end));
% sound(partirapoint,fe)



%Question 7:
% parole=[riennesertde,courir,ilfaut,partirapoint];
% sound(parole,fe)


