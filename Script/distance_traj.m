% DISTANCE TRAJ
% 
% Funzione per il calcolo della distanza media fra la retta della traiettoria del
% fireball e le rette delle osservazioni fatte da tutte le stazioni.
%
% Algoritmo tratto da: Borovicka, Bull. Astron. Inst. Czechosl 41 (1990),
% 391-396
%
% INPUT:
% Nstaz = numero delle stazioni osservative
% XM, YM, ZM = coordinate rettangolari geocentriche del punto iniziale della traiettoria (m)
% csi_R, eta_R, zeta_R = versori della traiettoria del bolide (adim.)
% X, Y, Z = coordinate delle stazioni osservative
% csi, eta, zeta = versori della traiettoria del bolide osservata dalle
% diverse stazioni al suolo.
% Numero = numero totale di posizioni osservate del bolide
%
% OUTPUT:
% D = distanza media fra la retta della traiettoria del
% fireball e le rette delle osservazioni fatte da tutte le stazioni (m)
% 
% Albino Carbognani, Versione del 15 gennaio 2019



function [D]=distance_traj(Nstaz, XM, YM, ZM, csi_R, eta_R, zeta_R, X, Y, Z, csi, eta, zeta, Numero)

D=0; d=0;

VM=[XM YM ZM];               % Vettore delle coordinate del punto iniziale della traiettoria del fireball

VersM=[csi_R eta_R zeta_R];  % Versore della traiettoria del fireball

for i=1:Nstaz

Vstaz=[X(i) Y(i) Z(i)];      % Vettore delle coordinate della stazione i-esima

    for j=1:length(csi{i})

          VersObs=[csi{i}(j) eta{i}(j) zeta{i}(j)]; % Versore della retta stazione i-esima, osservazione j-esima

          numeratore=abs(dot(VM-Vstaz, cross(VersM, VersObs)));

          denominatore=norm(cross(VersM, VersObs));
    
          d = d + numeratore/denominatore;

    end

D=D+d;

end

D=D/Numero;

end
