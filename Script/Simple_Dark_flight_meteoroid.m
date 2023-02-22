% DARK FLIGHT - Simplified Meteoroid impact point computation
%
% Modello semplificato della fase di dark flight di un meteoroide per ottenere una stima del punto di caduta al suolo. 
% L'algoritmo è tratto da: Z. Ceplecha, Bull. Astron. Inst. Czechosl. 38 (1987), 222-234.
%
% IPOTESI FISICHE DI BASE: La rotazione terrestre è trascurata, la velocità del vento è considerata nulla 
% e il modello di atmosfera considerato è quello esponenziale che meglio fitta il modello COESA del 1976. 
% Il moto del meteoroide avviene in un piano. Lx e Vx sono la distanza e la velocità orizzontale, mentre Ly e Vy sono le analoghe
% quantità verticali.
%
% DATI DI INPUT:
% data_path = stringa con il path del file da cui leggere i dati di ingresso 
% working_path = stringa con il path di salvataggio dei risultati
% fireball_name = stringa con il nome del fireball
%
% DATI LETTI DA FILE:
% lat_phi = Latitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radianti)
% long_lambda = Longitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radianti)
% a_R = Azimut (N -> E) del radiante apparente del meteoroide nell'ultimo punto osservato (in gradi ma trasformato in radianti)
% beta = Azimut (N -> E) della direzione verso cui si muove il bolide (in radianti)
% Zr = Distanza zenitale della traiettoria nell'ultimo punto osservato (in gradi ma trasformato in radianti)
% vT = Modulo della velocità del meteoroide nell'ultimo punto osservato in km/s         
% hT = Quota dell'ultimo punto osservato in km
% at = Accelerazione del meteoroide nell'ultimo punto osservato in km/s^2
% G_asintotico = Valore del coefficiente di resistenza aerodinamica gamma per numeri di Mach superiori a 4 (valore asintotico)
% Dfin = Rapporto M/S nel punto finale dal modello dinamico del meteoroide(kg/km^2);
% GS_su_M_fin = Prodotto Drag*S/M nel punto finale dalla cinematica del modello dinamico del meteoroide(km^2/kg);
%
% OUTPUT:
% I risultati vengono salvati nel file "Simple_Dark_flight_'fireball_name'.txt". 
% Il valore finale di Lx fornisce la distanza orizzontale (in km), fra la proiezione al suolo dell'ultimo
% punto osservato della traiettoria del fireball e il punto di caduta del meteoroide lungo la traiettoria percorsa. 
% Nel file di output vengono salvati anche i valori di latitudine e longitudine del punto di impatto al suolo.
%
% A. Carbognani - OAVdA
% Versione: 23 marzo 2019

function [] = Simple_Dark_flight_meteoroid(data_path, working_path, fireball_name)

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%    Simplified DARK FLIGHT - Meteoroid impact point computation      %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')


% Costanti
C=57.29577951;                              % Conversione da gradi a radianti e viceversa

% Inizializzazione dei dati sul fireball

disp('   ')
disp(' INITIALIZATION OF THE FIREBALL DATA  ')
disp('   ')

cd(data_path)
s0=[fireball_name '_Dynamics.txt'];
Dat(:,1)=load(s0, '-ascii');

lat_phi=Dat(1, 1)/C;                          % Latitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radiant
long_lambda=Dat(2, 1)/C;                      % Longitudine dell'ultimo punto osservato della traiettoria (in gradi ma trasformato in radianti)
a_R=Dat(3, 1)/C;                              % Azimut (N -> E) del radiante apparente del meteoroide nell'ultimo punto osservato (in gradi ma trasformato in radianti)
beta=a_R-pi;                                  % Azimut (N -> E) della direzione verso cui si muove il bolide (in radianti)
Zr=Dat(4, 1)/C;                               % Distanza zenitale della traiettoria nell'ultimo punto osservato (in gradi ma trasformato in radianti)
vT=Dat(5, 1);                                 % Modulo della velocità del meteoroide nell'ultimo punto osservato in km/s         
hT=Dat(6, 1);                                 % Quota dell'ultimo punto osservato in km
at=Dat(7, 1);                                 % Accelerazione del meteoroide nell'ultimo punto osservato in km/s^2
G_asintotico=Dat(8, 1);                       % Valore del coefficiente di resistenza aerodinamica gamma per numeri di Mach superiori a 4 (valore asintotico)
Dfin=Dat(9, 1)*10^6;                          % Rapporto M/S nel punto finale dal modello dinamico del meteoroide(kg/km^2);
GS_su_M_fin=Dat(10, 1)/10^6;                  % Prodotto Drag*S/M nel punto finale dalla cinematica del modello dinamico del meteoroide(km^2/kg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Integrazione numerica delle equazioni diff. per velocità e distanza   %                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' NUMERICAL INTEGRATION OF THE MOTION EQUATIONS  ')
disp('   ')

vx=vT*sin(Zr);                     % Velocità orizzontale iniziale lungo la traiettoria in km/s
vy=-vT*cos(Zr);                    % Velocità verticale iniziale in km/s
t_interval=(0:0.5:500);            % Vettore dei tempi di integrazione (in s)

% Calcolo dell'altezza di scala e della densità al suolo efficaci dell'atmosfera 
% in modo da fittare al meglio il modello COESA del 1976     

Height1=(0:0.1:hT); % Vettore delle quote da hT al suolo (km)

% Best fit di COESA con il modello atmosferico esponenziale

[rho0, Hs] = best_fit_COESA_model_dark_flight(Height1);

% Densità dell'aria (kg/km^3) alla quota iniziale osservata del dark flight 
% usando il modello esponenziale che meglio fitta il modello COESA del 1976.

rho_start=rho0*exp(1).^(-hT/Hs); 

% Velocità orizzontale km/s, velocità verticale km/s, distanza orizzontale km, distanza verticale km, densità aria alla quota iniziale kg/km^3;
init_cond = [vx, vy, 0, hT, rho_start]'; 

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
disp(' SAVE DARK FLIGHT DATA IN A .TXT FILE  ')
disp('   ')

fid = fopen([working_path '\Simple_Dark_flight_' fireball_name '.txt'],'w'); 
fprintf(fid, [' %% Simple dark flight model for the fireball: ' fireball_name '\n']);
fprintf(fid, ' %% Hypothesis: no Earth rotation; no winds; exponential atmosphere COESA 1976 model. \n');
fprintf(fid, ' %% Input data from atmospheric dynamical model. \n');
fprintf(fid, '    \n');
fprintf(fid, ' INPUT DATA \n');
fprintf(fid, '    \n');
fprintf(fid, ' %% Drag coefficient (adimensional):                    %10.2f \n', G_asintotico);
fprintf(fid, ' %% M/S for dark flight (kg/m^2):                       %10.2f \n', Dfin/10^6);
fprintf(fid, ' %% Drag*S/M for dark flight (m^2/kg):                  %10.6f \n', GS_su_M_fin*10^6);
fprintf(fid, ' %% Start height dark flight (km):                      %10.2f \n', hT);
fprintf(fid, ' %% Start velocity dark flight (km/s):                  %10.2f \n', vT);
fprintf(fid, ' %% Start acceleration dark flight (km/s^2):            %10.2f \n', at);
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

string1=['Simplified dark flight - Meteoroid height vs. horizontal distance_' fireball_name '.bmp'];
string1a=[fireball_name ' - Simplified dark flight - Meteoroid height vs. horizontal distance'];
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
disp(' END SIMPLE DARK FLIGHT COMPUTATION  ')
disp('   ')

end
