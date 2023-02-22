% Funzione vettoriale per il calcolo delle coordinate equatoriali di un corpo
% celeste a partire da quelle azimutali. L'azimut e' contato da
% nord verso est. Il tempo siderale calcolato dal programma e' quello
% medio, questo significa che l'ascensione retta degli astri al meridiano
% è quella media, quindi le coordinate equatoriali sono riferite all'equinozio 
% medio della data. 
%
% Input:
% LAT: latitudine del luogo in radianti
% AZ, H: altezza e azimut in radianti. Attenzione: Azimut contato da nord verso est.
% TS: tempo siderale medio locale in radianti
%
% Output:
% AR, DEC: coordinate equatoriali equinozio alla data in radianti
 
function [AR, DEC] = Azimut2EQ(LAT, AZ, H, TS)

 % Calcolo delle coordinate equatoriali (equinozio medio alla data)

  AZ=AZ-pi; % Passaggio da Azimut contato da Nord verso Est ad Azimut contato da Sud verso Ovest
 
  DEC=sin(LAT)*sin(H)-cos(LAT)*cos(H).*cos(AZ);

  DEC=asin(DEC);

  A1=(cos(H).*sin(AZ))./cos(DEC);

  A2=(cos(LAT)*sin(H)+sin(LAT)*cos(H).*cos(AZ))./cos(DEC);

  AH=atan2(A1,A2);

  AR=TS-AH;

 % Riduzione di AR fra 0 e 2 pi greco

 if AR>=0

    n=floor(AR/(2*pi));
    
 else 

    n=ceil(AR/(2*pi));
  
 end
    
 AR=AR-n*2*pi;

 if AR<0

    AR=AR+2*pi;

 end
 