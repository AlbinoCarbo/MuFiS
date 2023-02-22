% Funzione per il bin dei vettori osservati di quota, velocità e tempo
%
% INPUT: 
% height = vettore delle quote osservate (km)
% velocity = vettore delle velocità osservate (km/s)
% time = tempi osservati (s)
% bin = ampiezza del bin (s)
%
% OUTPUT:
% H = vettore quote binnate (km)
% V = vettore velocità binnate (km/s)
% T = tempi binnati (s)
%
% Versione 21 dicembre 2018 by Albino Carbognani (OAVdA)

function [XN, YN, ZN, H, V, T]=binning_vector(xn, yn, zn, height, velocity, time, bin)

N=length(time);
M=floor((max(time)-min(time))/bin)+1; % Numero di bin
tempo_min = min(time);
tempo_max = tempo_min+bin;

for j=1:M
    
    XX(j)=0; YY(j)=0; ZZ(j)=0; HH(j)=0; VV(j)=0; TT(j)=0; counter(j)=0; % Inizializzazione dei vettori binnati
    
    for i=1:N
        if time(i) <= tempo_max && time(i) >= tempo_min
        
          XX(j)=XX(j)+xn(i); YY(j)=YY(j)+yn(i); ZZ(j)=ZZ(j)+zn(i); HH(j)=HH(j)+height(i); VV(j)=VV(j)+velocity(i); TT(j)=TT(j)+time(i); 
          counter(j)=counter(j)+1;
        end
    end
    
    tempo_max=tempo_max+bin; tempo_min=tempo_min+bin; % Nuovi limiti del bin
    
end

% Calcolo valori medi dei vettori binnati
XN=XX./counter; YN=YY./counter; ZN=ZZ./counter;
H=HH./counter; V=VV./counter; T=TT./counter;

end