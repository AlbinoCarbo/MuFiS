% SUNLONG
% Albino Carbognani (OAVdA)
% Versione del 27 giugno 2017
%
% funzione per il calcolo della longitudine del Sole riferita
% all'equinozio 2000. Parametro formale: data giuliana. Il valore
% in output della longitudine del Sole è espresso in radianti.

function [lo]=sunlong(jd)

gr=57.29578; % costante di conversione gradi <--> radianti

% calcolo del secolo giuliano */

 t=(jd-2415020.0)/(36525.0); 

% calcolo della longitudime media del Sole */

 lo=(279.69668+36000.76892*t+0.0003025*t*t)/gr;

% calcolo della anomalia media del Sole */

 am=(358.47583+35999.04975*t-0.00015*t*t-0.0000033*t*t*t)/gr;

% calcolo equazione del centro del Sole */

 aux1=((1.91946-0.004789*t-0.000014*t*t)*sin(am))/gr;
 
 aux2=((0.020094-0.0001*t)*sin(2*am))/gr;
 
 aux3=(0.000293*sin(3*am))/gr;
 
 anno=1900+100*t;
 
 lo=lo+aux1+aux2+aux3-0.000243752*(anno-2000);

% riduzione della longitudine fra 0 e 2*pi */

 if(lo>=0)

    n=floor(lo/(2*pi));

 else 

    n=ceil(lo/(2*pi));
    
 end

 lo=lo-n*2*pi;

 if(lo<0)

    lo=lo+2*pi;
 end

