% DARK FLIGHT BASIC - Simplified Meteoroid impact point basic computation
%
% Modello semplificato della fase di dark flight di un meteoroide per ottenere una stima del punto di caduta al suolo senza passare per il
% computo del modello dinamico del meteoroide ma con i soli dati di triangolazione. 
% L'algoritmo è tratto da: Z. Ceplecha, Bull. Astron. Inst. Czechosl. 38 (1987), 222-234.
%
% IPOTESI FISICHE DI BASE: La rotazione terrestre è trascurata, la velocità del vento è considerata nulla 
% e il modello di atmosfera adottato è quello esponenziale con altezza di scala efficace pari a 7.64 km 
% (Il valore per l'altezza di scala si ottiene considerando che la densità dell'aria 
% diminuisce da 1200 g/m^3 a livello del mare a 0.125 g/m^3 a 70 km, che corrisponde a un 
% fattore 9600, indicando una altezza di scala "efficace" di 70/ln(9600) = 7,64 km).
% Il moto del meteoroide avviene in un piano. Lx e Vx sono la distanza e la velocità orizzontale, mentre Ly e Vy sono le analoghe
% quantità verticali. 
%
% NOTA: Questa funzione viene usata solo se si utilizza MuFiS in modalità "basic", ossia quando viene 
% usata la sola cinematica per caratterizzare il punto finale del bolide. Per questo motivo
% i risultati forniti sul dark flight e le coordinate del punto di caduta non sono sufficientemente 
% precisi per la ricerca di meteoriti al suolo. I valori sono solo indicativi!!
%
% DATI DI INPUT: 
% lat_phi, long_lambda: latitudine e longitudine del punto finale osservato della traiettoria (°)
% a_R = azimut (N -> E), del radiante apparente del bolide nell'ultimo punto osservato (°).
% Zr = distanza zenitale della traiettoria Zr nell'ultimo punto osservato.
% vT = modulo della velocità nell'ultimo punto osservato della traiettoria (km/s).
% hT = quota finale della traiettoria osservata (km). 
% aT = accelerazione nel punto finale osservato (km/s^2)
%
% OUTPUT:
% I risultati vengono salvati nel file "Simple_Dark_flight_'fireball_name'.txt". 
% Il valore finale di Lx fornisce la distanza orizzontale (in km), fra la proiezione al suolo dell'ultimo
% punto osservato della traiettoria del fireball e il punto di caduta del meteoroide lungo la traiettoria percorsa. 
% Nel file di output vengono salvati anche i valori di latitudine e longitudine del punto di impatto al suolo.
%
% A. Carbognani - OAVdA, versione del 21 febbraio 2019

function [] = Simple_Dark_flight_meteoroid_basic(working_path, fireball_name, lat_phi, long_lambda, a_R, Zr, vT, hT, aT)

disp('                                                                             ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%  Simplified DARK FLIGHT Basic - Meteoroid impact point basic computation  %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                             ')


% Costanti
C=57.29577951;                              % Conversione da gradi a radianti e viceversa
Hs=7.64;                                    % Altezza di scala efficace dell'atmosfera in km
rho_zero=1.225*10^9;                        % Densità dell'atmosfera al suolo in kg/km^3

% Inizializzazione dei dati sul fireball

disp('   ')
disp(' INITIALIZATION OF THE FIREBALL DATA  ')
disp('   ')

lat_phi=lat_phi/C;                            % Latitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radiant
long_lambda=long_lambda/C;                    % Longitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radianti)
a_R=a_R/C;                                    % Azimut (N -> E) del radiante apparente del meteoroide nell'ultimo punto osservato (in gradi ma trasformato in radianti)
beta=a_R-pi;                                  % Azimut (N -> E) della direzione verso cui si muove il bolide (in radianti)
Zr=Zr/C;                                      % Distanza zenitale della traiettoria nell'ultimo punto osservato (in gradi ma trasformato in radianti)
G_asintotico=0.58;                            % Valore del coefficiente di resistenza aerodinamica gamma per numeri di Mach superiori a 4 (valore asintotico)
rhof=rho_zero*exp(1).^(-hT/Hs);               % Densità atmosferica alla quota finale in kg/km^3
GS_su_M_fin=-aT/(rhof*vT*vT);                 % Prodotto Drag*S/M nel punto finale dalla cinematica del modello dinamico del meteoroide(km^2/kg);
Dfin=G_asintotico/GS_su_M_fin;                % Rapporto M/S nel punto finale dal modello dinamico del meteoroide(kg/km^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Integrazione numerica delle equazioni diff. per velocità e distanza   %                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' NUMERICAL INTEGRATION OF THE MOTION EQUATIONS  ')
disp('   ')

vx=vT*sin(Zr);                     % Velocità orizzontale iniziale lungo la traiettoria in km/s
vy=-vT*cos(Zr);                    % Velocità verticale iniziale in km/s
t_interval=(0:0.5:500);            % Vettore dei tempi di integrazione (in s)

% Velocità orizzontale km/s, velocità verticale km/s, distanza orizzontale km, distanza verticale km, densità aria alla quota iniziale kg/km^3;
init_cond = [vx, vy, 0, hT, rho_zero*exp(1).^(-hT/Hs)]'; 

[t, y] = ode45(@(t,Y) odefcn(t, Y, GS_su_M_fin, Hs), t_interval, init_cond);

v=sqrt(y(:,1).*y(:,1)+y(:,2).*y(:,2)); % Modulo della velocità di caduta (km/s)

Lx=y(:,3); % Vettore delle distanze orizzontali (km)
Ly=y(:,4); % Vettore delle distanze verticali (km)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dell'indice massimo con ancora la quota positiva %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=0;
for i=1:length(Ly)
    if Ly(i) >= 0
       N=N+1;
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate geografiche del punto di impatto %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_long=Lx(N)*sin(beta);   % Variazione di longitudine dal punto finale osservato (in km)
delta_lat=Lx(N)*cos(beta);    % Variazione di latitudine dal punto finale osservato (in km)

% Lunghezza del del grado di longitudine in km

grado_longitudine=2*pi*(6371*cos(lat_phi))/360;

% Lunghezza del grado di latitudine in km

grado_latitudine=2*pi*(6371)/360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apertura del file per il salvataggio dei risultati dell'integrazione numerica %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' SAVE DARK FLIGHT BASIC DATA IN A .TXT FILE  ')
disp('   ')

fid = fopen([working_path '\Simple_Dark_flight_basic_' fireball_name '.txt'],'w'); 
fprintf(fid, [' %% Simple dark flight basic model for the fireball: ' fireball_name '\n']);
fprintf(fid, ' %% Hypothesis: no Earth rotation; no winds; exponential atmosphere with H0 = 7.64 km. \n');
fprintf(fid, ' %% Input data from triangulation and Ceplecha cinematical model. \n');
fprintf(fid, '    \n');
fprintf(fid, ' INPUT DATA \n');
fprintf(fid, '    \n');
fprintf(fid, ' %% Fixed drag coefficient (adimensional):              %10.2f \n', G_asintotico);
fprintf(fid, ' %% M/S for dark flight (kg/m^2):                       %10.2f \n', Dfin/10^6);
fprintf(fid, ' %% Drag*S/M for dark flight (m^2/kg):                  %10.6f \n', GS_su_M_fin*10^6);
fprintf(fid, ' %% Start height dark flight (km):                      %10.2f \n', hT);
fprintf(fid, ' %% Start velocity dark flight (km/s):                  %10.2f \n', vT);
fprintf(fid, ' %% Start acceleration dark flight (km/s^2):            %10.2f \n', aT);
fprintf(fid, '    \n');
fprintf(fid, ' RESULTS \n');
fprintf(fid, '    \n');
fprintf(fid, ' %% Horizontal length dark flight (km):                 %10.2f \n', max(y(N, 3))); 
fprintf(fid, ' %% Flight time (s):                                    %10.2f \n', t_interval(N)); 
fprintf(fid, ' %% Geographical Latitude strewn field (degree):        %10.4f  \n', lat_phi*C+delta_lat/grado_latitudine);
fprintf(fid, ' %% Geographical Longitude strewn field (degree):       %10.4f \n', long_lambda*C+delta_long/grado_longitudine);
fprintf(fid, '    \n');
fprintf(fid, '%% quota (km)     Vx (km/s)       Vy (km/s)      V (km/s)        Lx (km)     tempo di caduta (s)\n');
fprintf(fid, '    \n');

for i=1:N
    fprintf(fid,'%10.3f \t %10.3f \t %10.3f \t %10.3f \t\t %10.2f \t %10.2f \n', y(i,4), y(i,1), y(i,2), v(i), y(i,3), t_interval(i));
end

% Chiusura del file di output
status = fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot della traiettoria della caduta %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' SAVE DARK FLIGHT TRAJECTORY PLOT IN .BMP FORMAT ')
disp('   ')

% Plot dello spostamento in avanti

string1=['Simplified dark flight basic - Meteoroid height vs. horizontal distance_' fireball_name '.bmp'];
string1a=[fireball_name ' - Simplified dark flight basic - Meteoroid height vs. horizontal distance'];
figure
plot(Lx(1:N), Ly(1:N), 'k.')
hold on
grid
xlabel('Horizontal distance (km)','FontSize',14)
ylabel('Altitude (km)','Fontsize',14)
title(string1a,'FontSize',14)
hold off

cd(working_path)
saveas(gcf, string1, 'bmp')   % Salva le figure in bmp

disp('   ')
disp(' END SIMPLE DARK FLIGHT BASIC COMPUTATION  ')
disp('   ')

end
