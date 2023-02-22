 % TEMPOSIDERALE
 % 
 % Albino Carbognani (OAVdA)
 % Versione del 14 giugno 2017
 %
 % Algoritmo: Di Lucio Angeletti, Pietro Giannone, Esercizi e complementi
 % di astronomia, Edizioni Nuova Cultura, Roma 2009. Pag. 21-22.
 %
 % Funzione per il calcolo del tempo siderale medio.
 % Il valore di TS restituito è espresso in radianti. 
 %
 % Parametri formali: 
 % jd=data giuliana per le 0h UT 
 % ut=ora di tempo universale per cui si vuole il ts,
 % lambda=longitudine del luogo in radianti (positiva verso est). 

 function [ts] = temposiderale(jd, ut, lambda)

 k1=1.002737908; % conversione ore solari --> ore siderali 

 k2=0.261799387; % conversione ore --> radianti 

 % calcolo del secolo giuliano 

 t=(jd-2451545.0)/(36525.0); 

 % tempo siderale medio a Greenwich per 0h tu (in radianti) 

 ts=k2*(6.697374558+(2400.051337)*t+(0.00002586222)*t*t);

 % tempo siderale medio alla longitudine lambda al tempo ut (in radianti)

 ts=ts+lambda+k1*k2*ut;

 % riduzione del tempo siderale fra 0 e 2*pi

 if ts>=0

    n=floor(ts/(2*pi));
    
 else 

    n=ceil(ts/(2*pi));
  
 end
    
 ts=ts-n*2*pi;

 if ts<0

    ts=ts+2*pi;

 end
 
 




