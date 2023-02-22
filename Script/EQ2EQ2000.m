% Funzione vettoriale per la trasformazione delle coordinate equatoriali 
% da equinozio alla data a J2000.0.
% Algoritmo: J. Meeus, Astronomia con il computer, Hoepli 1990, Cap.14
%
% Input: AR, DEC: coordinate alla data in radianti
% JDin: giurno giuliano dell'epoca alla data
%
% Output: AR2000, DEC2000: coordinate equatoriali al J2000.0 in radianti

function [AR2000, DEC2000] = EQ2EQ2000(AR, DEC, JDin)

CC=(180/pi)*3600; % Conversione secondi - d'arco radianti

JDfin=2451545;    % Giorno giuliano per l'epoca di riferimento J2000.0 (1 gennaio 2000 ore 12 UT)

tau0=(JDin-2415020.313)/36524.2199;

tau=(JDfin-JDin)/36524.2199;

% Costanti della precessione in secondi d'arco

a=(2304.250+1.396*tau0)*tau+0.302*tau^2+0.018*tau^3;
b=a+0.791*tau^2+0.001*tau^3;
c=(2004.682-0.853*tau0)*tau-0.426*tau^2-0.042*tau^3;

% Conversione da secondi d'arco a radianti

a=a/CC;
b=b/CC;
c=c/CC;

% Calcolo nuove coordinate J2000.0

A=cos(DEC).*sin(AR+a);
B=cos(c)*cos(DEC).*cos(AR+a)-sin(c)*sin(DEC);
C=sin(c)*cos(DEC).*cos(AR+a)+cos(c)*sin(DEC);

AR2000=b+atan2(A, B);

if AR2000 < 0
    AR2000=AR2000+2*pi;
end
    
    DEC2000=asin(C);


