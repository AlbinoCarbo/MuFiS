% CEPLECHA TWO STATIONS - Fireball trajectory computation for two stations. 
%
% L'algoritmo è tratto da: Z. Ceplecha, Bull. Astron. Inst. Czechosl. 38 (1987), 222-234.
% 
% Questa funzione è derivata dallo script "Trajectory_two_stations.m" e permette di calcolare il punto iniziale e finale
% della traiettoria osservata del fireball dalla prime due stazioni dell'elenco di tutte le stazioni che hanno osservato il fireball. 
%
% Le coordinate in output sono rettangolari geocentriche e servono per
% costruire la traiettoria media multi-stazione con il metodo di
% Borovicka, Bull. Astron. Inst. Czechosl. 41 (1990), 391-396.
%
% INPUT
% station_A = '1'; station_B = '2'; % Numeri d'ordine identificativi delle due stazioni
% R, raggio geocentrico della stazione (m); 
% h, quota della stazione (m); 
% Lat_geo, latitudine geocentrica della stazione (radianti); 
% Long, longitudine geocentrica della stazione (radianti)
% TS, tempo siderale della stazione (radianti); 
% Tempo, tempo delle osservazioni del fireball visto dalla stazione A (JD); 
% Dec, declinazione alla data delle posizioni del fireball (radianti); 
% AR,ascensione retta alla data delle posizioni del fireball (radianti) 
%
% OUTPUT
% Coordinate rettangolari geocentriche del punto iniziale e finale della
% traiettoria osservata dalla stazione A (in m):
% X_inA, Y_inA, Z_inA, X_finA, Y_finA, Z_finA
%
% FUNZIONI ESTERNE RICHIESTE:
% piano_osservazioni.m, funzione per il calcolo del piano che passa dalla
% stazione osservativa e contiene i versori di tutti i punti della
% traiettoria osservata del fireball.
%
% zero2pi.m, funzione che riduce gli angoli fra 0 e 2*pi.
%  
% Albino Carbognani - OAVdA
% Versione 06 Sept 2018

function [X_inA, Y_inA, Z_inA, X_finA, Y_finA, Z_finA] = Ceplecha_two_stations(working_path, fireball_name, station_A, station_B, R_A, h_A, Lat_geoA, Long_A, TS_A, Tempo_A, Dec_A, AR_A, R_B, h_B, Lat_geoB, TS_B, Dec_B, AR_B)

C=57.29577951; % Costante di conversione gradi <--> radianti

% Trasformazione del tempo da JD a secondi di tempo

Tempo_A=86400*(Tempo_A-min(Tempo_A));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del piano delle osservazioni per la stazione A %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XA, YA, ZA, csi_A, eta_A, zeta_A, a_A, b_A, c_A, d_A] = piano_osservazioni(R_A, h_A, Lat_geoA, TS_A, Dec_A, AR_A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del piano delle osservazioni per la stazione B %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XB, YB, ZB, csi_B, eta_B, zeta_B, a_B, b_B, c_B, d_B] = piano_osservazioni(R_B, h_B, Lat_geoB, TS_B, Dec_B, AR_B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate rettangolari geocentriche dei punti della traiettoria media per la stazione A (in metri) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_NA=eta_A*c_A-zeta_A*b_A;
b_NA=zeta_A*a_A-csi_A*c_A;
c_NA=csi_A*b_A-eta_A*a_A;
d_NA=-a_NA.*XA-b_NA.*YA-c_NA.*ZA;

n_osservazioni_A=length(AR_A);

% Intersezione di tre piani per ottenere le coordinate rettangolari goeocentriche dei punti osservati lungo la traiettoria
% media
for i=1:n_osservazioni_A
    A=[a_A b_A c_A; a_B b_B c_B; a_NA(i) b_NA(i) c_NA(i)];
    B=-[d_A; d_B; d_NA(i)];
    P_NA=A\B;
    X_NA(i)=P_NA(1);
    Y_NA(i)=P_NA(2);
    Z_NA(i)=P_NA(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate sferiche e delle quote per i punti della   %
% traiettoria media osservata dalla stazione A.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_NA=atan2(Y_NA, X_NA);                     % Tempo siderale dei punti osservati da A (in radianti), compreso fra -pi e +pi

if theta_NA <0
        theta_NA=2*pi+theta_NA;                 % Conversione a tempo siderale positivo (se necessario) 
end
    
lambda_NA=Long_A+(theta_NA-TS_A);               % Longitudine dei punti osservati da A (in radianti)
lambda_NA=zero2pi(lambda_NA);                   % Riduzione di lambda fra 0 e 2*pi
R_hA=sqrt(X_NA.^2+Y_NA.^2+Z_NA.^2); %R_hA=sqrt((X_NA./cos(theta_NA)).^2+Z_NA.^2);    % Distanza geocentrica dei punti osservati da A (raggio terrestre + quota, in metri)
phi_NA=asin(Z_NA./R_hA);                        % Latitudine geocentrica dei punti osservati da A (in radianti), fra -pi/2 e +pi/2
H_A=R_hA-geocradius(C*phi_NA, 'WGS84');         % Quota del bolide sulla superficie terrestre in metri

phi_NA_geog=(geoc2geod(C*phi_NA, geocradius(C*phi_NA, 'WGS84'), 'WGS84'))/C; % Latitudine geografica (o geodetica) dei punti osservati da A (in radianti)
phi_NA_geog_verticale=phi_NA_geog+(H_A.*(phi_NA-phi_NA_geog))./R_hA;         % Latitudine geografica dei punti osservati da A proiettati verticalmente sulla superficie (in radianti)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     TEST SULLA TRAIETTORIA OSSERVATA    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test per la riproduzuine della traiettoria
% Best_traiettoria calcolata - traiettoria osservata

theta_N=atan2(Y_NA-YA, X_NA-XA);              % Tempo siderale dei punti visti da A (in radianti), compreso fra -pi e +pi

if theta_N <0
        theta_N=2*pi+theta_N;                 % Conversione a tempo siderale positivo (se necessario) 
end

D_AR=AR_A-theta_N';
distanza=sqrt((X_NA-XA).^2+(Y_NA-YA).^2+(Z_NA-ZA).^2);   % Vettore delle distanze fra osservatore e fireball
D_Dec=Dec_A-(asin((Z_NA-ZA)./distanza))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Curva delle velocità per la stazione A.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vettore delle lunghezze percorse dal bolide fra la prima immagine e le
% successive (m)

for i=1:n_osservazioni_A
l_A(i)=sqrt((X_NA(1)-X_NA(i))^2+(Y_NA(1)-Y_NA(i))^2+(Z_NA(1)-Z_NA(i))^2);
end

% Calcolo della curva di velocità per la stazione A (m/s) 
Tempo_A1=Tempo_A';
v_A=gradient(l_A, Tempo_A1);                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Punti medi iniziale e finale della traiettoria osservata per la stazione A %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dati per il punto medio iniziale della traiettoria per la stazione A

X_inA=X_NA(1); Y_inA=Y_NA(1); Z_inA=Z_NA(1);

% Dati per il punto medio finale della traiettoria per la stazione A

X_finA=X_NA(n_osservazioni_A); Y_finA=Y_NA(n_osservazioni_A); Z_finA=Z_NA(n_osservazioni_A);

%==========================================================================

% Salvataggio dei dati numerici della traiettoria e velocità osservata dalla
% stazione A

fid0 = fopen([working_path '\' fireball_name '_Ceplecha_Trajectory_' station_A '_with_' station_B '.txt'],'w'); % Path di salvataggio del file di dati
fprintf(fid0, ' %%POSITIONAL DATA                                                                                              \n');
fprintf(fid0, ' %%Time (s)     Xgeoctrc (km)     Ygeoctrc (km)     Zgeoctrc (km)     lat (deg)      lon (deg)      height (km)           V (km/s)      Delta AR (deg)       Delta Dec (deg)     distanza osservatore (km)\n');
for i=1:n_osservazioni_A
fprintf(fid0,'%10.4f \t %10.4f \t\t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t\t %10.4f \t %10.4f \t\t %10.4f \t\t %10.2f\n', Tempo_A(i), X_NA(i)/1000, Y_NA(i)/1000, Z_NA(i)/1000, phi_NA_geog_verticale(i)*C, lambda_NA(i)*C, H_A(i)/1000, v_A(i)/1000, C*D_AR(i), C*D_Dec(i), distanza(i)/1000);
end

% Chiusura del file di output

status = fclose(fid0);

% Apertura dei file per il salvataggio dei risultati delle due stazioni A e B in file kml visualizzabili con Google Earth

filename1 = [working_path '\' fireball_name '_Ceplecha_Trajectory_' station_A '_with_' station_B '.kml'];
kmlwriteline(filename1, C*phi_NA_geog_verticale, C*lambda_NA, 'Color','yellow','Width', 10);

end

