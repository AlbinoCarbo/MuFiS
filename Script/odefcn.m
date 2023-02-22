% Funzione per la dinamica di un meteoroide in atmosfera nella fase di volo
% buio. Integra le equazioni del moto della resistenza aerodinamica.
%
% NOTA: il modello atmosferico per la densità dell'aria è quello esponenziale 
% che meglio fitta il modello COESA del 1976.
%
% Versione del 21 febbraio 2019 by Albino Carbognani

function [dYdt] = odefcn(t, Y, GD, Hs)

% Variabili:
         %Y(1); Velocità orizzontale Vx in km/s
         %Y(2); Velocità verticale Vy in km/s
         %Y(3); spostamento orizzontale in km
         %Y(4); quota in km
         %Y(5); densità dell'aria in kg/km^3
         
% Parametri: GD = Gamma*S/M_fin, Hs (altezza di scala dell'atmosfera)

g=9.81/1000; % Accelerazione di gravità in km/s^2

dYdt = [          
         -(GD)*Y(5).*sqrt(Y(1).^2+Y(2).^2).*Y(1);             % Equazione differenziale per Vx in funzione del tempo
         -(GD)*Y(5).*sqrt(Y(1).^2+Y(2).^2).*Y(2)-g;           % Equazione differenziale per Vy in funzione del tempo
         Y(1);                                                % Equazione differenziale per la distanza orizzontale in funzione del tempo 
         Y(2);                                                % Equazione differenziale per la quota in funzione del tempo
         -Y(2).*Y(5)/Hs;                                      % Equazione differenziale per la densità dell'aria in funzione del tempo
                    ];                                                        

end