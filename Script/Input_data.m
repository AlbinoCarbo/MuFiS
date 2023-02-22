% Funzione per l'input dei dati osservati nella funzione di triangolazione multistazione
% Può usare sia il nuovo formato adottato dal 1 aprile 2019, sia il vecchio formato, valido fino a IT20190222.
%
% INPUT:
% data_name = nome del file della stazione n-esima da cui leggere i dati
% input_format = tipo di input, 'P' (PRISMA) = nuovo formato, 'S' (standard) =
% vecchio formato.
% fireball_name = nome del fireball
%
% OUTPUT:
% julian_date = data giuliana degli istanti osservati del bolide
% Azi = azimut della traiettoria osservata in radianti (equinozio alla data)
% Alt = altezza sull'orizzonte della traiettoria osservata in radianti (equinozio alla data)
% Mag = magnitudine apparente osservata
%
% By A. Carbognani, 10 maggio (OAVdA)

function [julian_date, Azi, Alt, Mag]=Input_data(data_name, input_format, fireball_name)

% Costanti ausiliarie

C=57.29577951;   % Conversione da gradi a radianti e viceversa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nuovo formato di input adottato dal 1 aprile 2019 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = strcmp(input_format, 'P'); % Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
if tf == 1

fid = fopen(data_name,'rt');
indata = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',10); % Le prime 10 linee non vengono lette perché fanno parte dell'header
date=indata{1};    
julian_date=indata{2};
az=indata{3};
s_az=indata{4};
zd=indata{5};
s_zd=indata{6};
ra=indata{7};
s_ra=indata{8};
dec=indata{9};
s_dec=indata{10};
mag=indata{11};
s_mag=indata{12};
fclose(fid);

Azi = az/C;         % Azimut della traiettoria osservata dalla i-esima stazione (in radianti) 

% Bolidi per cui, nel formato di input PRISMA, la colonna n.5 è la distanza
% zenitale invece dell'altezza sull'orizzonte:

tf1 = strcmp(fireball_name, 'IT20190327');
tf2 = strcmp(fireball_name, 'IT20190405');
tf3 = strcmp(fireball_name, 'IT20190412');

if (tf1==1) || (tf2==1) || (tf3==1)

Alt = (90-zd)/C;   % Altezza della traiettoria osservata dalla i-esima stazione (in radianti) 

else
    
Alt = (zd)/C;      % Altezza della traiettoria osservata dalla i-esima stazione (in radianti)

end

Mag=mag;           % Magnitudine apparente osservata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vecchio formato di input valido fino al 31 marzo 2019 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else

  B(:,:)=load(data_name, '-ascii');
  julian_date=B(:,1);
  Azi = B(:, 4)/C;                   % Azimut della traiettoria osservata dalla i-esima stazione (in radianti)
  Alt = B(:, 5)/C;                   % Altezza della traiettoria osservata dalla i-esima stazione (in radianti)
  
  N=length(julian_date);
  Mag=zeros(1, N);                   % Inizializzazione vettore delle mag apparenti fittizie

  for i=1:N
  Mag(i)=99; % Magnitudine fittizia (il vecchio formato non ha la colonna delle magnitudini apparenti da leggere)
  end
  
end

end


