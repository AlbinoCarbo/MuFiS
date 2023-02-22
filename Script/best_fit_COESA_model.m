% Funzione per il calcolo della densità al suolo e dell'altezza di scala efficace 
% di un modello esponenziale dell'atmosfera in modo da fittare al meglio il modello COESA del 1976.
%
% FUNZIONI USATE:
% rho_atmosphere_COESA.m, Calcola la densità dell'atmosfera in funzione della quota 
% secondo il "Committee on Extension to the Standard Atmosphere" (COESA) model del 1976. 
% A questo scopo, al suo interno, usa la funzione "atmoscoesa" di Matlab. 
%
% INPUT:
% H = vettore delle quote (km)
%
% OUTPUT:
% Hs_best = altezza di scala efficace (km)
% rho_zero_best = densità aria alla quota 0 km (kg/km^3)
%
% By A. Carbognani (OAVdA), 21 febbraio 2019

function [rho_zero_best, Hs_best ] = best_fit_COESA_model(H)

% Costanti globali
Hs=7.64;                               % Altezza di scala efficace dell'atmosfera standard esponenziale (km)
rho_zero=1.225*10^9;                   % Densità dell'atmosfera standard esponenziale al suolo (kg/km^3)

% Soppressione del Warning: The Jacobian at the solution is ill-conditioned, and some model parameters may not be estimated well (they are not identifiable).  Use caution in making predictions. 
warning('off', 'stats:nlinfit:IllConditionedJacobian');

Model_atmosphere=@(q, x)(q(1)*exp(1).^(-x/q(2))); % Modello atmosferico esponenziale

density_COESA=rho_atmosphere_COESA(H); % Densità dell'atmosfera nel modello COESA (kg/km^3)

% Valori iniziali stimati dei parametri p(1) e p(2) 
startingVals = [rho_zero Hs];  

% Fit non lineare dei parametri
[coefEsts, R1, J1, CovB1, MSE1] = nlinfit(H, density_COESA, Model_atmosphere, startingVals);

rho_zero_best=coefEsts(1); % Densità efficace atmosfera alla quota 0

Hs_best=coefEsts(2); % Altezza di scala efficace dell'atmosfera (km)

fit_density=Model_atmosphere(coefEsts, H); % kg/km^3

% Plot per il test del fit della densità atmosferica
figure
plot(H, rho_atmosphere_COESA(H), 'r.', 'MarkerSize', 16)    % Valori COESA della densità atmosferica
hold on
plot(H, fit_density, 'k-', 'LineWidth', 3)       % Valori della densità atmosferica con il modello esponenziale corretto
hold on
plot(H, Model_atmosphere(startingVals, H), 'b-', 'LineWidth', 3)       % Valori della densità atmosferica con il modello esponenziale originale
grid
legend('COESA model','Modified exponential model', 'Original exponential model')
ylabel('Air density (kg/km^3)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title('Test plot: atmospheric models - Exponential vs COESA (1976)','FontSize',20)
hold off

end

