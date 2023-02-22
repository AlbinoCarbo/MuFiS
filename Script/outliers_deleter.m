% Outliers deleter
%
% Funzione per la cancellazione degli elementi di un vettore che sono a più
% di 3 sigma dalla mediana
%
% Albino Carbognani (OAS)
% Versione del 6 maggio 2020

function [clean]=outliers_deleter(vector)

sigma=std(vector); % Deviazione standard del vettore
m=median(vector);  % Mediana del vettore

j=1;
N=length(vector);

for i=1:N
    if abs(vector(i)-m) < 3*sigma
        vettore(j)=vector(i);
        j=j+1;
    end
end

clean=vettore; % Nuovo vettore ripulito dagli elementi oltre i 3 sigma dalla mediana

end

    