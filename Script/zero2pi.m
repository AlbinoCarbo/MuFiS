% ZERO2PI
% Albino CArbognani (OAVdA)
% Versione del 2 settembre 2018
%
% Funzione per la riduzione di un angolo angolo fra 0 e 2*pi
% Input: valore di un angolo in radianti
% Output: angolo ridotto fra 0 e 2*pi

function [y] = zero2pi(x)

 if x >= 0

    n = floor(x/(2*pi)); % La funzione floor arrotonda un numero al numero intero più piccolo successivo
                       % Esempi: floor(3.5) = 3; floor(-7/2) = -4;
 else 

    n = ceil(x/(2*pi));  % La funzione ceil arrotonda un numero al numero intero più grande successivo
                       % Esempi: floor(3.5) = 4; floor(-7/2) = -3;
 end
    
 y = x-n*2*pi;

 if y < 0

    y = y+2*pi;

 end