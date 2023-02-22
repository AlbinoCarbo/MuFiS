% Funzione rho_atmosphere_COESA
%
% Calcola la densità dell'atmosfera in funzione della quota 
% secondo il "Committee on Extension to the Standard Atmosphere" (COESA)
% model del 1976. A questo scopo usa la funzione "atmoscoesa" di Matlab.
%
% NOTA sull'input/output della funzione atmoscoesa
% [T, a, P, Rho] = atmoscoesa(h)
%
% h = array delle quote in m 
% T = array of temperatures in kelvin.
% a = array speeds of sound, in m/s. The function calculates speed of sound using a perfect gas relationship.
% P = array of pressures, in Pascal.
% Rho = array of densities, in kg/m^3. The function calculates density using a perfect gas relationship.
%
% INPUT:
% heights = array delle quote (km)
%
% OUTPUT:
% rho = array della densità dell'atmosfera in kg/km^3
%
% By A. Carbognani, 15 febbraio 2019


function [ rho ] = rho_atmosphere_COESA(heights)

% Soppressione del Warning: A height is outside the range of 0 meters and 84852 meters. Output values will be extrapolated for those heights.
warning('off', 'aero:atmoscoesa:tooLowWarn');

% Calcolo della densità dell'atmosfera alla quota di input secondo il modello COESA (1976), in kg/m^3
% Bisogna moltiplicare la quota di input della funzione "rho_atmosphere_COESA" per 1000 così viene trasformata in metri.
[T, a, P, rho] = atmoscoesa(heights*1000); 

% Conversione della densità atmosferica in kg/km^3
rho=(10^9)*rho;                               

end

