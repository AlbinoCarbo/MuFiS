% PIANO OSSERVAZIONI
% Versione del 21 giugno 2017
% Albino Carbognani (OAVdA)
%
% Funzione per il calcolo del piano delle osservazioni che passa dalla
% stazione osservativa
%
% Input:
% R = raggio terrestre geocentrico della stazione in m
% h = quota della stazione in m
% Lat_geo = latitudine geocentrica in radianti
% TS = tempo siderale della osservazione del bolide in radianti
% Dec = declinazioni della traiettoria del bolide in radianti (equinozio alla
% data)
% AR = ascensione retta della traiettoria del bolide in radianti (equinozio
% alla data)
%
% Output:
% X, Y, Z = Coordinate rettangolari geocentriche della stazione in m
% a, b, c, d = costanti del piano passante per la stazione e contenente i
% versori delle osservazioni della traiettoria (adimensionali)
% csi, eta, zeta = versori geocentrici della traiettoria osservata dalla
% stazione (in radianti)


function [X, Y, Z, csi, eta, zeta, a, b, c, d] = piano_osservazioni(R, h, Lat_geo, TS, Dec, AR)

% Coordinate rettangolari geocentriche della stazione
X=(R+h)*cos(Lat_geo)*cos(TS);
Y=(R+h)*cos(Lat_geo)*sin(TS);
Z=(R+h)*sin(Lat_geo);

% Versori geocentrici della traiettoria osservata
csi=cos(Dec).*cos(AR);
eta=cos(Dec).*sin(AR);
zeta=sin(Dec);

a1=sum(csi.*eta)*sum(eta.*zeta)-sum(eta.*eta)*sum(csi.*zeta);
b1=sum(csi.*eta)*sum(csi.*zeta)-sum(csi.*csi)*sum(eta.*zeta);
c1=sum(csi.*csi)*sum(eta.*eta)-(sum(csi.*eta))^2;
d1=sqrt(a1^2+b1^2+c1^2);

% Componenti del vettore unitario perpendicolare al piano medio delle
% osservazioni della stazione
a=a1/d1;
b=b1/d1;
c=c1/d1;

% Distanza fra il piano passante per la stazione (contenente le osservazioni), e il centro della Terra

d=-(a*X+b*Y+c*Z); 
