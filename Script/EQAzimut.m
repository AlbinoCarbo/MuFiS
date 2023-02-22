% EQAzimut
%
% Albino Carbognani (OAVdA)
% Versione del 16 giugno 2017
%
% Funzione per il calcolo delle coordinate azimutali di un corpo
% celeste a partire da quelle equatoriali. L'azimut e' contato da
% nord verso est. Il tempo siderale e' quello medio, questo significa 
% che l'ascensione retta degli astri al meridiano
% è quella media, quindi le coordinate equatoriali da convertire in
% azimutali devono essere riferite all'equinozio medio della data.

% Input:
% LAT, latitudine del luogo in radianti
% AR, DEC, coordinate equatoriali alla data in radianti
% ts, tempo siderale locale in radianti
%
% Output:
% H, altezza sull'orizzonte in radianti
% AZ, azimut contato da nord verso est in radianti

function [H, AZ] = EQAzimut(LAT, AR, DEC, ts)

% Calcolo delle coordinate azimutali

  H=sin(LAT)*sin(DEC)+cos(LAT)*cos(DEC).*cos(ts-AR);

  H=asin(H);

  A1=(-cos(LAT)*sin(DEC)+sin(LAT)*cos(DEC).*cos(ts-AR))./cos(H);

  A2=(cos(DEC).*sin(ts-AR))./cos(H);

  AZ=atan2(A2,A1);   % Azimut compreso fra -pi e +pi radianti
  
  AZ=AZ+pi;           % Trasformazione in azimut contato da nord verso est (in radianti)