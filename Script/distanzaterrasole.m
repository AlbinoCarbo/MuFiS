% TERRASOLE
% Albino Carbognani (OAVdA)
% Versione del 27 giugno 2017
%
% Funzione per il calcolo della distanza Terra-Sole in UA. 
% Parametri formali: giorno giuliano.

function [R]=distanzaterrasole(jd)

gr=57.29578; % costante di conversione gradi <--> radianti

% calcolo della data in secoli giuliani

 t=(jd-2415020.0)/(36525.0);

% calcolo della anomalia media del Sole

 am=(358.47583+35999.04975*t-0.00015*t*t-0.0000033*t*t*t)/gr;

% calcolo della eccentricita' dell'orbita terrestre

 e=0.01675104-0.0000418*t-0.000000126*t*t;

% risoluzione iterativa della equazione di Keplero

 ae=am;

 for i=1:10
   
    ae=am+e*sin(ae);
    
 end
   
% calcolo della distanza Terra-Sole in UA 

 R=(1.0000002*(1-e*cos(ae)));