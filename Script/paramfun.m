% Funzione per l'integrazione delle equazioni differenziali del moto 
% che descrivono la dinamica di un meteoroide che cade in atmosfera (CONSTANT-ABLATION METHOD)
%
% IPOTESI FISICHE DI BASE: 1) Il meteoroide deve mantenere la sua forma durante
% l'ablazione 2) La superfice terrestre è considerata piatta. 
% Modello di Kalenichenko, A&A 448,(2007). 
% 
% NOTA: la densità dell'aria è espressa da una funzione esponenziale in modo che 
% fitti al meglio il modello COESA del 1976 (modello non isotermico).
%
% Versione del 21 febbraio 2019 by Albino Carbognani

function pos = paramfun(x, tspan)

G=x(1);            % Gamma, coefficiente di drag
D0=x(2);           % M/S, rapporto massa - sezione meteoroide all'infinito, kg/km^2
s=x(3);            % sigma, coefficiente di abalzione (s/km)^2
Vinf=x(4);         % Velocità all'infinito del meteoroide km/s
z=x(5);            % Coseno dell'inclinazione media della traiettoria
Hs=x(6);           % Altezza di scala efficace dell'atmosfera scelta in modo che fitti il modello COESA del 1976.
init_cond=x(7:9);  % Condizioni iniziali per velocità, densità aria e altezza


% Solution, Kalenichenko, A&A 448,(2007)
f =@(t, Y) [          
         -(G/D0)*Y(2).*Y(1).*Y(1).*exp(1)^(-(s/6)*(Y(1).*Y(1)-Vinf*Vinf)); % Equazione differenziale dV/dt=-(G/D0)*Rho*V^2*exp(-sigma/6*(V^2-Vinf^2)) per la velocità V in funzione del tempo
         Y(2).*Y(1)*z/Hs;                                                  % Equazione differenziale dRho/dt=dRho*V*cos(z)/Hs per la densità dell'aria Rho in funzione del tempo
         -Y(1)*z];                                                         % Equazione differenziale per la quota H in funzione del tempo

     [~,pos] = ode45(f,tspan,init_cond);                                   % Runge-Kutta 4th/5th order ODE solver

end