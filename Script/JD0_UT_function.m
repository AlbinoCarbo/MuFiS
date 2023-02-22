% Funzione temporale per la pipeline dei fireball
%
% FUNZIONI ESTERNE RICHIESTE:
%
% Input_data.m, funzione per l'input dei dati dalla i-esima stazione
% osservativa
%
% INPUT: 
% data_path = path dei dati 
% fireball_name = nome del fireball
% input_format = tipo di input, 'P' (PRISMA) = nuovo formato, 'S' (standard) =
% vecchio formato.
%
% OUTPUT: 
% jd0 = UT medio del fireball 
% UT = giorno giuliano dell'evento alle 0 UT
%
% Albino Carbognani (OAVdA)
%
% Last update: 10 maggio 2019

function [jd0, UT] = JD0_UT_function(data_path, fireball_name, input_format)

cd(data_path)    % Legge nella cartella dei dati osservati

s1=[fireball_name '_1.txt'];                 % Legge il file dati della prima stazione osservativa (che deve essere il migliore per estensione temporale)

[JD, Azi, Alt, mag]=Input_data(s1, input_format, fireball_name); % Estrazione JD osservati dalla prima stazione osservativa

% Giorno giuliano medio della osservazione del bolide
JD_medio=mean(JD);

% Estrazione parte decimale del JD medio e calcolo dell'ora media del bolide in UT
integ=floor(JD_medio);
UT=12.0+(JD_medio-integ)*24.0; % Ora media del fireball in UT

if UT > 24 % UT è > 24 h se il fireball è osservato al mattino, se sì questo ciclo if-end normalizza l'ora.
    UT=UT-24;
end

jd0=JD_medio-(UT/24); % JD per le 0 UT

end