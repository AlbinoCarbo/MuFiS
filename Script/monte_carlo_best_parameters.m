% Funzione per la stima della incertezza dei parametri di best fit del modello 
% dinamico dei bolidi con il metodo Monte Carlo.
%
% Input:
% Monte_Carlo_Number= numero di cicli Monte Carlo da calcolare
%
% init_cond = matrice delle condizioni iniziali (i parametri D0, s e Vinf
% vengono fatti variare casualmente al massimo del +/- 3*Frac, dove Frac è 
% la frazione di variazione media dello scenario di velocità e quota).
% 
% Tempo1=vettore temporale dei dati osservati
% soln=matrice dei valori osservati di velocità, densità aria e quote
% (vengono fatti variare casualmente al max del +/- 3*DH; 3*DV).
% lb, ub=constrain sui parametri di best fit
% options=opzioni per lsqcurvefit
% DV=scarto medio valori velocità (km/s)
% DH=scarto medio valori quote (km)
%
% Output:
% Sigma_G=deviazione standard della media coefficiente di drag
% sigma_D0=deviazione standard della media rapporto massa - sezione meteoroide all'infinito (kg/m^2)
% sigma_s=deviazione standard della media coefficiente di ablazione (s^2/km^2)
% sigma_Vinf=deviazione standard della media velocità all'infinito del meteoroide (km/s)
% sigma_V_start=deviazione standard della media velocità iniziale del meteoroide (km/s) 
%
% Albino Carbognani
% Versione del 22 aprile 2020

function [sigma_G, sigma_D0, sigma_s, sigma_Vinf, sigma_V_start] = monte_carlo_best_parameters(Monte_Carlo_Number, init_cond, Tempo1, soln, lb, ub, options, DV, DH)

% Inizializzazione dei vettori dei parametri
G=zeros(1, Monte_Carlo_Number);          % Gamma, coefficiente di drag
D0=zeros(1, Monte_Carlo_Number);         % M/S, rapporto massa - sezione meteoroide all'infinito, kg/m^2
s=zeros(1, Monte_Carlo_Number);          % sigma, coefficiente di ablazione (s^2/km^2)
Vinf=zeros(1, Monte_Carlo_Number);       % Velocità all'infinito del meteoroide km/s 
V_start=zeros(1, Monte_Carlo_Number);    % Velocità iniziale del meteoroide in km/s 

% Creazione matrici degli scenari Monte Carlo per velocità e quote
% N.B. randn è un numero casuale con distribuzione normale standard, ossia con valore medio
% zero e deviazione standard 1.

Matrix_DV=DV*randn(length(soln(:,1)), Monte_Carlo_Number); 
Matrix_DH=DH*randn(length(soln(:,3)), Monte_Carlo_Number); 

FracV=DV/mean(soln(:,1)); % Frazione di variazione dello scenario delle velocità
FracH=DH/mean(soln(:,3)); % Frazione di variazione dello scenario delle quote
Frac=(FracV+FracH);       % Frazione di variazione media dello scenario osservato

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inizio ciclo Monte Carlo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Monte_Carlo_Number

% Valori randomizzati dei vettori osservati di velocità, densità dell'aria e quota
soln1(:,1)=soln(:,1)+Matrix_DV(:, i);       % Velocità randomizzata del fireball in funzione del tempo
soln1(:,2)=soln(:,2);                       % Densità dell'aria in funzione della quota (non viene randomizzata)
soln1(:,3)=soln(:,3)+Matrix_DH(:, i);       % Altezza randomizzata del fireball in funzione del tempo
    
% Valori randomizzati delle condizioni iniziali cui è consentito variare max +/- 3*Frac
init_cond(2)=init_cond(2)+Frac*init_cond(2)*randn;              % D0=M/S, rapporto massa - sezione meteoroide all'infinito, kg/km^2
init_cond(3)=init_cond(3)+Frac*init_cond(3)*randn;              % s, coefficiente di abalzione (s/km)^2
init_cond(4)=init_cond(4)+DV*randn;                             % Vinf, velocità all'infinito del meteoroide km/s

% Calcolo parametri di best fit per lo scenario Monte Carlo
[pbest, ~, ~, ~, ~, ~, ~,] = lsqcurvefit(@paramfun, init_cond, Tempo1, soln1, lb, ub, options);

% Inizializzazione vettori dei parametri 
G(i)=pbest(1);                   
D0(i)=pbest(2)/(10^6);           
s(i)=pbest(3);                   
Vinf(i)=pbest(4);                 
V_start(i)=pbest(7);              

end

% Calcolo delle deviazioni standard della media dei parametri di best fit
sigma_G=std(G)/sqrt(Monte_Carlo_Number);
sigma_D0=std(D0)/sqrt(Monte_Carlo_Number);
sigma_s=std(s)/sqrt(Monte_Carlo_Number);
sigma_Vinf=std(Vinf)/sqrt(Monte_Carlo_Number);
sigma_V_start=std(V_start)/sqrt(Monte_Carlo_Number);

end
