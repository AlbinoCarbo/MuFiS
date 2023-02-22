% TRAJECTORY MULTI STATIONS - Fireball trajectory computation for N stations. 
%
% L'algoritmo per la triangolazione della traiettoria di un fireball da N stazioni contemporaneamente 
% è tratto da: J. Borovicka, Bull. Astron. Inst. Czechosl. 41 (1990), 391-396.
%
% FUNZIONI ESTERNE RICHIESTE:
%
% Input_data.m, funzione per l'input dei dati dalla i-esima stazione
% osservativa
%
% temposiderale.m, funzione per il calcolo del tempo siderale medio locale in radianti al momento della
% osservazione del fireball.
% 
% Azimut2EQ.m, funzione per il passagio da coordinate azimutali (azimut contato da nord verso
% est) a coordinate equatoriali alla data.
%
% EQ2EQ2000.m, funzione per il passaggio dalla coordinate equatoriali all'equinozio medio della data all'equinozio J2000.0
%
% zero2pi.m, funzione per la riduzione di un angolo fra 0 e 2*pi.
%
% Ceplecha_two_stations.m, funzione derivata da "Trajectory_two_stations.m". Permette di calcolare una traiettoria preliminare del fireball 
% usando i dati delle prime due stazioni dell'elenco da caricare. Per questo motivo le prime due stazioni devono essere quelle con le osservazioni migliori.  
% 
% ATTENZIONE: Se il numero di stazioni osservative è superiore a due, la traiettoria preliminare viene successivamente raffinata minimizzando la somma delle distanze fra le rette delle direzioni 
% osservate e la retta che descrive la traiettoria del fireball. In caso contrario si mantiene la soluzione della triangolazione di Ceplecha.
%
% inter_line_plane.m, funzione per il calcolo delle coordinate del punto di intersezione fra una retta e un piano.
%
% distance_traj.m, funzione per il calcolo della somma delle distanze fra la retta della traiettoria del
% fireball e le rette delle osservazioni fatte da tutte le stazioni.
%
% plane_generator.m, funzione per il calcolo delle costanti di un piano a partire da due vettori giacenti sul piano e un punto
% da cui il piano deve passare.
%
% piano_osservazioni.m, funzione usata all'interno di "Ceplecha_two_stations.m" per il calcolo del piano che passa dalla
% stazione osservativa e contiene i versori di tutti i punti della traiettoria osservata del fireball.
%
% EQAzimut.m, funzione che trasforma da coordinate equatoriali alla data in azimutali con azimut contato da nord verso est.
%
% Ceplecha_model_velocity.m, Funzione per la stima preliminare della velocità all'infinito (con incertezza), del meteoroide con il   
% modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h). Utile per la stima preliminare dell'orbita con i soli dati di
% triangolazione, senza passare dal modello dinamico che, se le osservazioni non sono accurate, potrebbe non convergere.
%
% INPUT
% data_path = path dei dati osservati, ossia "fireball_name\Data"
% working_path = path della cartella per l'output dei risultati della triangolazione, ossia "fireball_name\Trajectory"
% fireball_name = nome del fireball
% jd0 = giorno giuliano alle 0 UT dell'osservazione del fireball
% UT = ora in UT dell'osservazione del fireball
%
% OUTPUT
% Azimut_fin = Coordinate azimutali del radiante apparente geocentrico del fireball (radianti)
% fireball_total_duration = durata totale del fireball (s)
% Nstaz = numero di stazioni che hanno osservato il fireball
% AR_R = ascensione retta alla data del radiante apparente geocentrico del fireball (radianti)
% Dec_R = declinazione alla data del radiante apparente geocentrico del fireball (radianti)
% Lat_fin = latitudine del punto finale osservato del fireball (gradi) 
% Long_fin = longitudine del punto finale osservato del fireball (gradi)
% radiant_uncertainty = incertezza sulla posizione in cielo del radiante (gradi)
% V_infty = guess della velocità all'infinito nel modello di Ceplecha (km/s)
% DV_infty = deviazione standard della velocità all'infinito nel modello di Ceplecha (km/s)
% quota_finale = altezza finale del bolide (km)
% d_quota = incertezza media della quota (km)
% V_finale = velocità finale modello di Ceplecha in km/s
% Azimut_finale della traiettoria del bolide in gradi da N -> E
% inclinazione_finale della traiettoria del bolide in gradi
% Afin = accelerazione nel punto finale della traiettoria stimata usando il
% modello di Ceplecha delle velocità (km/s^2).
%
% SALVATAGGIO SU FILE 
% 1-file di testo con tutti i dati di triangolazione della traiettoria, 
% 2-Tutte le figure visualizzate a schermo in formato .bmp 
% 3-Tre file .kml visualizzabili con Google Earth che riportano la traiettoria proiettata al suolo e in quota. 
%
% Nota: i file .kml con la traiettoria in quota sono due. Il primo visualizza i punti osservati della traiettoria media. 
% Il secondo visualizza anche le linee che dai punti della traiettoria arrivano al suolo.
%
% Albino Carbognani - OAVdA
% Versione 12 maggio 2020

function [Azimut_fin, fireball_total_duration, Nstaz, AR_R, Dec_R, Lat_fin, Long_fin, radiant_uncertainty, V_infty, DV_infty, quota_finale, d_quota, V_finale, Azimut_finale, inclinazione_finale, Afin] = Trajectory_multi_stations(data_path, working_path, fireball_name, input_format, jd0, UT)

% clear all

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%        BOROVICKA TRAJECTORY - Path computation for N stations       %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')

% Costanti ausiliarie

C=57.29577951;   % Conversione da gradi a radianti e viceversa

% DATI DI INPUT PER LE STAZIONI OSSERVATIVE E DATI OSSERVATI DELLA TRAIETTORIA DEL FIREBALL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lettura del file con i dati delle stazioni che hanno osservato il fireball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_path)                                       % Legge nella cartella dei dati osservati

s0=['Stazioni_' fireball_name '.txt'];
A(:,:)=load(s0, '-ascii');

Lat = A(:,1)/C;                                     % Latitudine geografica (o geodetica) della stazione (radianti)
Long = A(:,2)/C;                                    % Longitudine geografica (o geodetica) della stazione (radianti)
h = A(:,3);                                         % Quota della stazione s.l.m. (m)
station = A(:,4);                                   % Numeri identificativi delle stazioni
Nstaz=length(station);                              % Numero di stazioni che hanno osservato il fireball
Lat_geo=(geod2geoc(C*Lat, h, 'WGS84'))/C;           % Latitudine geocentrica delle stazioni (radianti)
R=geocradius(C*Lat_geo, 'WGS84');                   % Raggio geocentrico della stazione (m)
TS=temposiderale(jd0, UT, Long);                    % Vettore del tempo siderale medio locale di osservazione del fireball nelle diverse stazioni (radianti)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load dei dati della traiettoria osservata dalle stazioni %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocazione dei vettori (velocizza i calcoli)
X=zeros(1, Nstaz); Y=zeros(1, Nstaz); Z=zeros(1, Nstaz);
a=zeros(1, Nstaz); b=zeros(1, Nstaz); c=zeros(1, Nstaz); d=zeros(1, Nstaz);

Nobs=0;                                 % Inizializzazione contatore del numero delle osservazioni fatte sul bolide
tempo_iniziale=zeros(1, Nstaz);         % Inizializzazione vettore dei tempi di osservazione iniziali per ciascuna stazione

for i=1:Nstaz % Lettura dei dati osservati dalla singole stazioni e inizializzazione delle matrici di tipo cell array.

s1=[fireball_name '_' num2str(station(i)) '.txt'];

[julian_date, Azi, Alt, mag]=Input_data(s1, input_format, fireball_name); % Funzione per l'input dei dati osservati dalla i-esima stazione


Mag_{i}=mag;                            % Matrice (cell array), formata da vettori colonna con lunghezza diversa, contenenti le mag apparenti osservate del fireball
Tempo_{i}=julian_date;                  % Matrice (cell array), formata da vettori colonna con lunghezza diversa, contenenti i tempi osservati del fireball (in giorni giuliani)
                                        % NB: questi vettori temporali servono per il calcolo della traiettoria iniziale provvisoria

tempo_iniziale(i)=min(julian_date);     % Vettore contenente i tempi iniziali di tutte le sequenze temporali osservate dalle stazioni (in giorni giuliani)
time{i}=julian_date;                    % Matrice (cell array), formata da vettori colonna con lunghezza diversa, contenenti i tempi osservati del fireball (in giorni giuliani)

Distanza_zenitale_{i}=pi/2-Alt;         % Vettore contenente le distanze zenitali del bolide osservato dalla i-esima stazione (in radianti). 
                                        % Questi dati sono usati per correggere le magnitudini apparenti dall'assorbimento atmosferico.

% Conversione da coordinate altazimutali (in radianti) ad equatoriali alla
% data media (in radianti)

[AR, Dec] = Azimut2EQ(Lat(i), Azi, Alt, TS(i));

Nobs=Nobs+length(AR);                  % Numero progressivo delle osservazioni compiute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% CREAZIONE MATRICI DI AR E DEC DELLE OSSERVAZIONI %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrici (cell array), formate da vettori colonna, con lunghezza diversa, contenenti AR e Dec 
% alla data media della traiettoria osservata dalla stazione (in radianti)
%
% Per avere la i-esima colonna della matrice AR: AR_{i} i=numero stazione

AR_{i}=AR;
Dec_{i}=Dec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calcolo del piano medio passante per le osservazioni e la stazione i-esima  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X0, Y0, Z0, csi0, eta0, zeta0, a0, b0, c0, d0] = piano_osservazioni(R(i), h(i), Lat_geo(i), TS(i), Dec, AR);

% Vettori con le coordinate rettangolari geocentriche X, Y, Z delle stazioni (in m)

X(i)=X0; Y(i)=Y0; Z(i)=Z0;           

% Inizializzazione dei vettori a, b, c, d aventi come componenti le costanti del piano medio passante per la stazione e contenente i
% versori delle osservazioni della traiettoria 

a(i)=a0; b(i)=b0; c(i)=c0; d(i)=d0;  

% csi0, eta0, zeta0 = versori geocentrici della traiettoria media osservata dalla
% stazione i-esima (in radianti)                                     
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% CREAZIONE MATRICI DEI VERSORI DELLE OSSERVAZIONI %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrici (cell array), formate da vettori colonna, con lunghezza diversa, contenenti i versori 
% geocentrici della traiettoria osservata dalla stazione (in radianti)
%
% Per avere la i-esima colonna della matrice csi: csi{i} i=numero stazione.
% Idem per le altre.

csi{i}=csi0;
eta{i}=eta0;
zeta{i}=zeta0;

clear julian_date Azi Alt mag AR Dec csi0 eta0 zeta0

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trasformazione della matrice del tempo delle osservazioni da JD a secondi (t=0 s è il primo punto osservato della traiettoria)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nstaz
   
    time{i}=86400*(time{i}-min(tempo_iniziale)); 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Triangolazione iniziale "alla Ceplecha" con le prime due stazioni dell'elenco che devono avere le migliori osservazioni disponibili.
%                        Se possibile le due stazioni dovranno vedere il bolide da due direzioni opposte.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' STARTING TRIANGULATION WITH TWO STATIONS -  CEPLECHA MODEL  ')
disp('   ')

station_A = '1'; station_B = '2'; % Numeri d'ordine identificativi delle prime due stazioni

[X_inA, Y_inA, Z_inA, X_finA, Y_finA, Z_finA] = Ceplecha_two_stations(working_path, fireball_name, station_A, station_B, R(1), h(1), Lat_geo(1), Long(1), TS(1), Tempo_{1}, Dec_{1}, AR_{1}, R(2), h(2), Lat_geo(2), TS(2), Dec_{2}, AR_{2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valori guess iniziali per i parametri della traiettoria rettilinea del fireball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinate guess del punto iniziale (in m)

XM_best=X_inA; 
YM_best=Y_inA; 
ZM_best=Z_inA; 

% Versori guess del vettore della traiettoria del fireball. 
% L'orientamento va dal punto finale verso il punto iniziale osservato (adimensionali)

% csi_R=(X_inA-X_finA)/abs(X_inA-X_finA); % Algoritmo originale di Borovicka
% eta_R=(Y_inA-Y_finA)/abs(X_inA-X_finA);
% zeta_R=(Z_inA-Z_finA)/abs(X_inA-X_finA);

vector_length=sqrt((X_inA-X_finA).^2+(Y_inA-Y_finA).^2+(Z_inA-Z_finA).^2);

csi_R_best=(X_inA-X_finA)/vector_length;       % Mio algoritmo
eta_R_best=(Y_inA-Y_finA)/vector_length;
zeta_R_best=(Z_inA-Z_finA)/vector_length;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Ottimizzazione dei parametri della traiettoria del fireball        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE BEST FIREBALL TRAJECTORY  ')
disp('   ')

% Distanza media per punto osservato fra la traiettoria guess del fireball e i versori di
% tutte le osservazioni (in m)

[Diniziale]=distance_traj(Nstaz, XM_best, YM_best, ZM_best, csi_R_best, eta_R_best, zeta_R_best, X, Y, Z, csi, eta, zeta, Nobs);

if Nstaz > 2 % WARNING: L'operazione di ottimizzazione della traiettoria viene fatta solo se il numero delle stazioni è maggiore di 2.
             % In caso contrario viene mantenuta la triangolazione di Ceplecha.

% Funzione da minimizzare variando i parametri della retta della
% traiettoria del fireball

Best_fit_traj=@(b)distance_traj(Nstaz, b(1), b(2), b(3), b(4), b(5), b(6), X, Y, Z, csi, eta, zeta, Nobs);

% Parametri iniziali della traiettoria del fireball

b_guess = [XM_best YM_best ZM_best csi_R_best eta_R_best zeta_R_best];

% Minimizzazione della somma delle distanze fra la retta del fireball e le
% rette di tutti i punti osservati

disp('   ')
disp(' OPTIMIZATION OF THE FIREBALL TRAJECTORY  ')
disp('   ')

% Plot per visualizzare l'operazione di minimizzazione della distanza:
options = optimset('PlotFcns',@optimplotfval); 

b_min = fminsearch(Best_fit_traj, b_guess, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametri della retta del fireball che minimizzano la distanza media per
% punto osservato.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Best punto della traiettoria del fireball

XM_best=b_min(1);  
YM_best=b_min(2);
ZM_best=b_min(3);

% Normalizzazione a 1 dei best coseni direttori della traiettoria del fireball

csi_R_best=b_min(4)/sqrt(b_min(4)^2+b_min(5)^2+b_min(6)^2);
eta_R_best=b_min(5)/sqrt(b_min(4)^2+b_min(5)^2+b_min(6)^2);
zeta_R_best=b_min(6)/sqrt(b_min(4)^2+b_min(5)^2+b_min(6)^2);

end

% Media delle distanze fra la best traiettoria del fireball e i versori di
% tutte le osservazioni (in m).
% NB: in un mondo ideale deve essere: Dfinale < Diniziale!

[Dfinale]=distance_traj(Nstaz, XM_best, YM_best, ZM_best, csi_R_best, eta_R_best, zeta_R_best, X, Y, Z, csi, eta, zeta, Nobs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate equatoriali del radiante apparente geocentrico medio del fireball %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTA: Il radiante apparente è il punto di intersezione fra la direzione di
% provenienza del fireball e la sfera celeste.

Dec_R=asin(zeta_R_best); % Declinazione geocentrica del fireball in radianti

XR=csi_R_best/cos(Dec_R);
YR=eta_R_best/cos(Dec_R);

AR_R=atan2(YR, XR);  % Ascensione retta geocentrica del fireball in radianti (compresa fra -pi e +pi)

if AR_R < 0
    AR_R=AR_R+2*pi;  % Riduzione dell'AR ad un valore sempre positivo compreso fra 0 e 2*pi
end

% Coordinate equatoriali J2000.0 del radiante apparente geocentrico del fireball calcolato con la traiettoria media

[AR2000_R, Dec2000_R] = EQ2EQ2000(AR_R, Dec_R, jd0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate rettangolari geocentriche dei punti della traiettoria media %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE RECTANGULAR GEOCENTRIC COORDINATES OF THE OBSERVED POINTS ALONG FIREBALL TRAJECTORY  ')
disp('   ')

a_best=zeros(1,Nstaz); b_best=zeros(1,Nstaz); c_best=zeros(1,Nstaz); d_best=zeros(1,Nstaz); % Pre-allocazione del vettore (velocizza i calcoli)

for i=1:Nstaz
    
    % Costanti del piano che passa per la i-esima stazione e la best traiettoria del
    % fireball
    
    [a_best(i), b_best(i), c_best(i), d_best(i)]=plane_generator(X(i)-XM_best, Y(i)-YM_best, Z(i)-ZM_best, csi_R_best, eta_R_best, zeta_R_best, XM_best, YM_best, ZM_best);
    
    ap=zeros(1,length(csi{i})); bp=zeros(1,length(csi{i})); cp=zeros(1,length(csi{i})); dp=zeros(1,length(csi{i})); % Pre-allocazione del vettore (velocizza i calcoli)
    xt=zeros(1,length(csi{i})); yt=zeros(1,length(csi{i})); zt=zeros(1,length(csi{i}));                             % Pre-allocazione del vettore (velocizza i calcoli)
    
    for j=1:length(csi{i})
             
    % Costanti del piano contenente la j-esima osservazione fatta dalla i-esima stazione
    
    [ap(j), bp(j), cp(j), dp(j)]=plane_generator(a_best(i), b_best(i), c_best(i), csi{i}(j), eta{i}(j), zeta{i}(j), X(i), Y(i), Z(i));
    
    % Coordinate rettangolari geocentriche dei punti della traiettoria
    % media per la i-esima stazione
     
    [xt(j), yt(j), zt(j)]=inter_line_plane(ap(j), bp(j), cp(j), dp(j), XM_best, YM_best, ZM_best, csi_R_best, eta_R_best, zeta_R_best);
    
    end
    
    % NOTA: Le coordinate Xgeoctrc (km), Ygeoctrc (km), Zgeoctrc (km) che si trovano nel file di output sono ottenute intersecando 
    % la retta della best traiettoria media con il piano contenente la j-esima osservazione della i-esima stazione e ortogonale al piano 
    % contenete la traiettoria media e la i-esima stazione. In pratica sono le coordinate del punto osservato che giace sulla traiettoria media. 
    % Infatti se si guarda la traiettoria nello spazio 3D è perfettamente rettilinea.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRICI DELLE COORDINATE DEI PUNTI OSSERVATI DELLA TRAIETTORIA DEL FIREBALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrici (cell array), formate da vettori colonna, con lunghezza diversa, contenenti le coordinate 
% geocentriche della traiettoria osservata dalle stazioni (in m)
%
% Per avere la i-esima colonna della matrice XT: XT{i} i=numero stazione
    
XT{i}=xt;
YT{i}=yt;
ZT{i}=zt;
    
clear ap bp cp dp xt yt zt

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo delle coordinate sferiche e delle quote per i punti della   %
% traiettoria media osservati dalla stazioni.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE HEIGHT, LATITUDE AND LONGITUDE OF THE OBSERVED POINTS ALONG FIREBALL TRAJECTORY  ')
disp('   ')

for i=1:Nstaz

    H=zeros(1, length(XT{i}));                     % Pre-allocazione del vettore delle quote osservate dalla i-esima stazione
    lambda=zeros(1, length(XT{i}));                % Pre-allocazione del vettore delle longitudini osservate dalla i-esima stazione
    phi_geog_verticale=zeros(1, length(XT{i}));    % Pre-allocazione del vettore delle latitudini osservate dalla i-esima stazione
    
    for j=1:length(XT{i})
    
       theta=atan2(YT{i}, XT{i});                  % Tempo siderale dei punti osservati dalla i-esima stazione (in radianti), compreso fra -pi e +pi
       
       [theta] = zero2pi(theta);                   % Conversione a tempo siderale compreso fra 0 e 2*pi 
           
       lambda=Long(i)+(theta-TS(i));               % Longitudine dei punti osservati dalla i-esima stazione (in radianti)

       [lambda]=zero2pi(lambda);                   % Riduzione di lambda fra 0 e 2*pi

       R_h=sqrt((XT{i}./cos(theta)).^2+ZT{i}.^2);  % Distanza geocentrica dei punti osservati dalla i-esima stazione (raggio terrestre + quota, in metri)

       phi=asin(ZT{i}./R_h);                       % Latitudine geocentrica dei punti osservati dalla i-esima stazione (in radianti), fra -pi/2 e +pi/2

       H=R_h-geocradius(C*phi, 'WGS84');           % Quota del bolide sulla superficie terrestre in metri

       phi_geog=(geoc2geod(C*phi, geocradius(C*phi, 'WGS84'), 'WGS84'))/C; % Latitudine geografica (o geodetica) dei punti osservati dalla i-esima stazione (in radianti)

       phi_geog_verticale=phi_geog+(H.*(phi-phi_geog))./R_h;              % Latitudine geografica dei punti osservati dalla i-esima stazione proiettati verticalmente sulla superficie (in radianti)

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRICI DELLE COORDINATE GEOGRAFICHE E DELLE QUOTE DEI PUNTI DELLA
% TRAIETTORIA DEL FIREBALL VISTA DALLE DIVERSE STAZIONI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrici (cell array), formate da vettori colonna, con lunghezza diversa, contenenti le coordinate 
% geografiche e le quote della traiettoria osservata dalle stazioni
%
% Per avere la i-esima colonna della matrice H: H{i} i=numero stazione
    
    HM{i}=H;
    lambdaM{i}=lambda;
    phi_geog_verticaleM{i}=phi_geog_verticale;
    
    clear H phi_geog_verticale lambda phi R_h theta
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo della matrice della magnitudine assoluta del fireball %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcola la magnitudine assoluta solo se il formato di input è quello di
% PRISMA

tf = strcmp(input_format, 'P'); % Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
if tf == 1

disp('   ')
disp(' COMPUTE FIREBALL ABSOLUTE MAG  ')
disp('   ')

for i=1:Nstaz
    
 % Calcolo coordinate rettangolari traiettoria bolide vista dalla stazione
 % i-esima
    
 % Lat(i) = Latitudine geografica (o geodetica) della i-esima stazione (radianti)
 % Long(i) = Longitudine geografica (o geodetica) della i-esima stazione (radianti)
 % h(i) = Quota della i-esima stazione s.l.m. (m)
 
[xNorth,yEast,zDown] = geodetic2ned(C*phi_geog_verticaleM{i},C*lambdaM{i},HM{i},C*Lat(i),C*Long(i),h(i),referenceEllipsoid);

% Vettore delle distanze fra la i-esima stazione e il fireball osservato
% (metri)

distance=sqrt(xNorth.^2+yEast.^2+zDown.^2)'; 

% Magnitudine assoluta del bolide visto dalla i-esima stazione

h_met=10^5;    % Quota di riferimento del bolide per il calcolo della magnitudine assoluta (100 km)

% Calcolo magnitudine assoluta del bolide. 
% NOTA: La mag apparente viene corretta per l'airmass assumendo un coefficiente di assorbimento in V di 0.1 mag/airmass.

mag_assoluta_bolide_{i}=(Mag_{i}-0.1*cos(Distanza_zenitale_{i}))-5*log10(distance/h_met);

clear distance

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Calcolo della matrice della velocità                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE FIREBALL VELOCITY  ')
disp('   ')

for i=1:Nstaz
    
    L=zeros(1, length(XT{i})); % Pre-allocazione del vettore delle lunghezze percorse
    
    for j=1:length(XT{i})
        % Vettore delle lunghezze percorse dal bolide fra la prima immagine e le
        % successive riprese dalla i-esima stazione(m)
        L(j)=sqrt((XT{i}(1)-XT{i}(j))^2+(YT{i}(1)-YT{i}(j))^2+(ZT{i}(1)-ZT{i}(j))^2);
    end
    
    % Calcolo vettore velocità per la stazione i-esima (m/s) 
    v=gradient(L, time{i}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRICE DELLE VELOCITA' OSSERVATE DEL FIREBALL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matrice (cell array), formata da vettori colonna, con lunghezza diversa, contenenti le
% velocità osservate del fireball dalle diverse stazioni.
%
% Per avere la i-esima colonna della matrice V: V{i} i=numero stazione
    
    V{i}=v;
        
    clear L v
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               OUTPUT DEI RISULTATI SULLA TRAIETTORIA                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' PLOT AND SAVE FIGURES IN .BMP FORMAT ')
disp('   ')

% ******************************************
% PLOT DELLA TRAIETTORIA MEDIA IN 3D
% ******************************************

% Plot della traiettoria geocentrica in 3D (in km)
string4=[fireball_name ' - 3D geocentric trajectory '];
figure
for i=1:Nstaz
plot3(XT{i}/1000, YT{i}/1000, ZT{i}/1000, 'k.', 'color', rand(1,3), 'MarkerSize', 16)
hold on
end
grid
xlabel('Geocentric X (km)','FontSize',20)
ylabel('Geocentric Y (km)','FontSize',20)
zlabel('Geocentric Z (km)','FontSize',20)
title(string4,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string4, 'bmp')   % Salva le figure in bmp

% ******************************************
% PLOT DELLA QUOTA VS TEMPO
% ******************************************

% Plot della quota (in km) vs tempo (in s)
string3=[fireball_name ' - Height vs time '];
figure
for i=1:Nstaz
plot(time{i}, HM{i}/1000, 'k.', 'color', rand(1,3), 'MarkerSize', 16, 'DisplayName', sprintf('Observed heights from station number %d', i))
legend('-DynamicLegend');
legend('show');
hold on
end
grid
%legend('Observed heights from mean trajectory')
xlabel('Time (s)','FontSize',20)
ylabel('Height (km)','FontSize',20)
title(string3,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string3, 'bmp')   % Salva le figure in bmp

% ******************************************
% PLOT DELLA VELOCITA' VS TEMPO
% ******************************************

% Plot della velocità (in km/s) vs tempo (in s)
string2=[fireball_name ' - Velocity vs time '];
figure
for i=1:Nstaz
plot(time{i}, V{i}/1000, 'k.', 'color', rand(1,3), 'MarkerSize', 16, 'DisplayName', sprintf('Observed velocity from station number %d', i))
legend('-DynamicLegend');
legend('show');
hold on
end
grid
xlabel('Time (s)','FontSize',20)
ylabel('Velocity (km/s)','FontSize',20)
title(string2,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string2, 'bmp')   % Salva le figure in bmp

% ********************************
% PLOT DELLA VELOCITA' VS QUOTA 
% ********************************

% Plot della quota (in km) vs velocità (in km/s)
string1=[fireball_name ' - Velocity vs Height '];
figure
for i=1:Nstaz
plot(HM{i}/1000, V{i}/1000, 'k.', 'color', rand(1,3), 'MarkerSize', 16, 'DisplayName', sprintf('Observed velocity from station number %d', i))       % Valori osservati della velocità
legend('-DynamicLegend');
legend('show');
hold on
end
grid
ylabel('Velocity (km/s)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title(string1,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string1, 'bmp')   % Salva le figure in bmp

% ******************************************************************
% PLOT DELLA MAG ASSOLUTA E APPARENTE (PER LA 1° STAZIONE) VS TEMPO
% ******************************************************************

% Plotta la mag assoluta e apparente vs tempo solo se il formato di input è quello di PRISMA
% Nel formato standard di input dei dati questo plot non c'è.

tf = strcmp(input_format, 'P'); % Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
if tf == 1

% Plot della mag assoluta vs tempo (in s)
string0=[fireball_name ' - Absolute mag vs time '];
figure
for i=1:Nstaz
plot(time{i}, mag_assoluta_bolide_{i}, 'k.', 'color', rand(1,3), 'MarkerSize', 16, 'DisplayName', sprintf('Observed absolute mag from mean trajectory %d', i))
set(gca,'Ydir','reverse')
legend('-DynamicLegend');
legend('show');
hold on
end
grid
xlabel('Time (s)','FontSize',20)
ylabel('Absolute Mag (mag)','FontSize',20)
title(string0,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string0, 'bmp')   % Salva le figure in bmp

% Plot della mag apparente vs tempo (in s) per la prima stazione
string=[fireball_name ' - Apparent mag vs time '];
figure
for i=1:Nstaz
plot(time{1}, Mag_{1}, 'k.', 'MarkerSize', 16)
set(gca,'Ydir','reverse')
hold on
end
grid
legend('Observed apparent mag from the first station')
xlabel('Time (s)','FontSize',20)
ylabel('Apparent Mag (mag)','FontSize',20)
title(string,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string, 'bmp')   % Salva le figure in bmp

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Salvataggio su file dei dati numerici della traiettoria e velocità osservata
% dalle diverse stazioni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' SAVE TRAJECTORY RESULTS IN .TXT FILES  ')
disp('   ')

fid0 = fopen([working_path '\' fireball_name '_Trajectory_' num2str(Nstaz) '_stations_raw.txt'],'w'); % Path di salvataggio del file di dati
fprintf(fid0, ['%% Fireball:' fireball_name '\n']);
fprintf(fid0, ' %%   \n');
fprintf(fid0, ' %%POSITIONAL DATA                                                                                              \n');
fprintf(fid0, ' %%Time (s)     Xgeoctrc (km)     Ygeoctrc (km)     Zgeoctrc (km)     lat (deg)      lon (deg)      height (km)           V (km/s)     N stazione\n');

for i=1:Nstaz
    for j=1:length(time{i})
        fprintf(fid0,'%10.4f \t %10.4f \t\t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t\t %10.4f \t\t %1.0f \n', time{i}(j), XT{i}(j)/1000, YT{i}(j)/1000, ZT{i}(j)/1000, C*phi_geog_verticaleM{i}(j), C*lambdaM{i}(j), HM{i}(j)/1000, V{i}(j)/1000, i);
    end
end

% Chiusura del file di output
status = fclose(fid0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ordinamento temporalmente crescente dei risultati  %
% NB: Legge, ordina e sovrascrive il file precedente %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1=[fireball_name '_Trajectory_' num2str(Nstaz) '_stations_raw.txt'];
D(:,:)=load(s1, '-ascii');

E=sortrows(D, 1); % Ordinamento della matrice D in senso crescente secondo la prima colonna (che corrisponde agli istanti di osservazione espressi in secondi)

% Calcolo del vettore dei coseni della distanza zenitale della best-traiettoria
% del fireball

LB=zeros(1, length(E(:, 1))-1); % Pre-allocazione dei vettori
HB=zeros(1, length(E(:, 1))-1);

for i=1:length(E(:, 1))-1
LB(i)=sqrt((E(1, 2)-E(i+1, 2))^2+(E(1, 3)-E(i+1, 3))^2+(E(1, 4)-E(i+1, 4))^2);
HB(i)=E(1, 7)- E(i+1, 7);
end

COS_Zr=HB./LB;

trajectory_length=LB(length(LB));                                    % Lunghezza della traiettoria percorsa in atmosfera dal bolide (km)
radiant_uncertainty=C*atan(((Dfinale)/1000)/trajectory_length);      % Stima della incertezza angolare del radiante del bolide (gradi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dell'inclinazione e dell'azimut del radiante geocentrico  %
% del bolide nel punto finale osservato                             %               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TS_fin=temposiderale(jd0, UT, E(length(E(:,1)), 6)/C);  % Tempo siderale medio locale nel punto finale osservato (radianti)

% Calcolo delle coordinate azimutali del radiante apparente geocentrico visto dal punto finale
% osservato del bolide.

[H_R, Az_R]=EQAzimut(E(length(E(:,1)), 5)/C, AR_R, Dec_R, TS_fin);

Inclinazione_fin=H_R; Azimut_fin=Az_R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costanti geometriche della traiettoria per la Terra sferica (Ceplecha, M&PS 40, ReVelle, 2005)
% Servono per calcolare l'inclinazione della traiettoria del bolide sulla superficie terrestre
% e la lunghezza del cammino osservato.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re=6378;                             % Raggio equatoriale della Terra in km
hmax=max(E(:,7));                    % Quota iniziale del fireball in km
AA=2*Re*COS_Zr(2)+2*hmax*COS_Zr(2);
BB=2*Re;  
CC=-2*Re*hmax-hmax*hmax;
incl_angle=abs(90-C*acos((AA/2-1)./(BB/2+E(:,7))));            % Inclinazione locale della traiettoria del bolide sulla superficie terrestre (gradi)
length_atmo=AA/2-sqrt((AA*AA)/4+CC+E(:,7).*E(:,7)+BB*E(:,7));  % Lunghezza del cammino osservato del bolide in atmosfera (km)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funzione per la stima preliminare della velocità all'infinito del meteoroide con il   
% modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h)  
% Utile per la stima preliminare dell'orbita con i soli dati di
% triangolazione, senza passare dal modello dinamico che, se le
% osservazioni non sono accurate, potrebbe non convergere.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V_infty, DV_infty, fit_velocity, Afin, tipo_modello]=Ceplecha_model_velocity(E(:, 7), E(:, 8), Inclinazione_fin);

% Plot del modello di Ceplecha o lineare, quota (in km) vs velocità (in km/s)
string4=[fireball_name ' - Velocity vs Height' tipo_modello ' model '];
figure
plot(E(:, 7), E(:, 8), 'r.', 'MarkerSize', 16)          % Valori osservati della velocità
hold on
plot(E(:, 7), fit_velocity, 'k-', 'LineWidth', 3)       % Curva delle velocità con il modello di Ceplecha (1961)
grid
legend('Observed values', [tipo_modello ' model'])
ylabel('Velocity (km/s)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title(string4,'FontSize',20)
hold off

cd(working_path)
saveas(gcf, string4, 'bmp')   % Salva le figure in bmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid1 = fopen([working_path '\' fireball_name '_Trajectory_' num2str(Nstaz) '_stations_ordered.txt'],'w'); % Path di salvataggio del file di dati

fprintf(fid1, [' %% Fireball: ' fireball_name '\n']);
fprintf(fid1, ' %% JD0: %10.6f \n', jd0);
fprintf(fid1, ' %% UT:  %10.6f \n', UT);
fprintf(fid1, ' %%    \n');
fprintf(fid1, ' %% Number of total observed points:                             %10.0f \n', length(E(:, 1)));
fprintf(fid1, ' %%    \n');
fprintf(fid1, ' %% Start height (km):                                             %10.1f \n', E(1,7));
fprintf(fid1, ' %% End height (km):                                               %10.1f \n', E(length(E(:, 7)),7));
fprintf(fid1, ' %% Atmospheric trajectory length (km):                            %10.1f \n', trajectory_length);
fprintf(fid1, ' %% Lat. and Long start (gradi):                                      %10.4f     %10.4f \n', E(1,5), E(1,6));
fprintf(fid1, ' %% Lat. and Long end (gradi):                                        %10.4f     %10.4f \n', E(length(E(:, 5)),5), E(length(E(:, 6)),6));
fprintf(fid1, ' %%    \n');
fprintf(fid1, ' %% Average trajectory inclination on the horizon (degree):           %10.2f \t pm  %10.2f\n', mean(incl_angle), std(incl_angle));
fprintf(fid1, ' %% Average meteoroid velocity (km/s):                               %10.1f \t\t pm %10.1f\n', mean(E(:, 8)), std(E(:, 8)));
fprintf(fid1, ' %% Guess pre-atmospheric velocity with the %s model (km/s):   %10.1f \t\t pm %10.1f\n', tipo_modello, V_infty, DV_infty);
fprintf(fid1, ' %% Guess final velocity with the %s model (km/s):             %10.1f \n', tipo_modello, fit_velocity(length(E(:,7))));
fprintf(fid1, ' %% Guess final acceleration with the %s model (km/s^2):       %10.1f \n', tipo_modello, Afin);
fprintf(fid1, ' %%    \n');
fprintf(fid1, ' %% Equatorial coordinates to date of the apparent geocentric radiant of the fireball computed with the average trajectory\n');
fprintf(fid1, ' %% RA (degree)     Dec (degree) \n');
fprintf(fid1, ' %% %10.4f \t %10.4f \n', AR_R*C, Dec_R*C);
fprintf(fid1, ' %%    \n');    
fprintf(fid1, ' %% Equatorial coordinates J2000.0 of the apparent geocentric radiant of the fireball computed with the average trajectory\n');
fprintf(fid1, ' %% RA (degree)     Dec (degree) \n');
fprintf(fid1, ' %% %10.4f \t %10.4f \n', AR2000_R*C, Dec2000_R*C);
fprintf(fid1, ' %%    \n');    
fprintf(fid1, ' %% Azimuthal coordinates of the apparent geocentric radiant of the fireball computed with the average trajectory from the last observed point\n');
fprintf(fid1, ' %% Az (degree)     Inclination (degree) \n');
fprintf(fid1, ' %% %10.2f \t %10.2f \n', Azimut_fin*C, Inclinazione_fin*C);
fprintf(fid1, ' %%    \n');    
fprintf(fid1, ' %% Estimate of spatial and angular uncertainties on the fireball trajectory: \n');
fprintf(fid1, ' %% Average distance between the observed points of the fireball trajectory and the trajectory computed with two-station  (km):  %10.4f \n', (Diniziale)/1000);  
fprintf(fid1, ' %% Average distance between the observed points of the fireball trajectory and the trajectory computed with N stations (km):    %10.4f \n', (Dfinale)/1000);
fprintf(fid1, ' %% Angular uncertainty on the radiant (degree):                                                                                 %10.4f \n', radiant_uncertainty);
fprintf(fid1, ' %%    \n');       
fprintf(fid1, ' %% POSITIONAL DATA                                                                                              \n');
fprintf(fid1, ' %%Time (s)     Xgeoctrc (km)     Ygeoctrc (km)     Zgeoctrc (km)     lat (deg)      lon (deg)      height (km)           V (km/s)       Inc. angle (degree)     Atmospheric path (km)      N stazion\n');

    for j=1:length(E(:,1))
        
        fprintf(fid1,'%10.7f \t %10.4f \t\t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t\t %10.4f \t\t %10.2f \t\t\t %10.2f \t\t\t\t %1.0f \n', E(j,1), E(j,2), E(j,3), E(j,4), E(j,5), E(j,6), E(j,7), E(j,8), incl_angle(j), length_atmo(j)-length_atmo(1), E(j,9));
        
    end

% Chiusura del file di output ordinato in senso temporalmente crescente
status = fclose(fid1);

fireball_total_duration=E(length(E(:,1)),1); % Durata totale in s del fireball

Lat_fin=E(length(E(:,1)), 5);  % Latitudine del punto finale osservato del fireball (in gradi) 

Long_fin=E(length(E(:,1)), 6); % Longitudine del punto finale osservato del fireball (in gradi) 

quota_finale=E(length(E(:, 7)),7); % Quota finale del bolide in km

d_quota=(Dfinale)/1000; % Incertezza media sulla quota (km)

V_finale =fit_velocity(length(E(:,7))); % Velocità finale del bolide nel modello di Ceplecha (km/s)

Azimut_finale = Azimut_fin*C; % Azimut della traiettoria nel punto finale (gradi) 

inclinazione_finale = incl_angle(length(E(:,7))); % Inclinazione della traiettoria nel punto finale (gradi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Salvataggio dei risultati della triangolazione in un file kml visualizzabile con Google Earth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' SAVE TRAJECTORY RESULTS IN .KML FILES FOR GOOGLE EARTH  ')
disp('   ')

% Syntax: kmlwriteline(filename,latitude,longitude)
% kmlwriteline(filename,latitude,longitude,altitude)

% ********************************************************************************
% Path delle icone che compariranno nel file kml visualizzabile con Google Earth 
% ********************************************************************************

% Icona locale (quando si riavvia Google Earth non viene caricata l'icona)
% iconFilename1 ='C:\Users\demo\Documents\MATLAB\Fireball\Multiple_Trajectories\Icone_Google_Earth\placemark_circle_highlight.png';
% Icona remota (quando si riavvia Google Earth viene caricata l'icona)
iconFilename1 = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png';

% Icona locale (quando si riavvia Google Earth non viene caricata l'icona)
% iconFilename2 = 'C:\Users\demo\Documents\MATLAB\Fireball\Multiple_Trajectories\Icone_Google_Earth\star.png';
% Icona remota (quando si riavvia Google Earth viene caricata l'icona)
iconFilename2 = 'http://maps.google.com/mapfiles/kml/shapes/star.png';

% Icona locale (quando si riavvia Google Earth non viene caricata l'icona)
% iconFilename3 = 'C:\Users\demo\Documents\MATLAB\Fireball\Multiple_Trajectories\Icone_Google_Earth\ylw-blank-lv.png';
% Icona remota (quando si riavvia Google Earth viene caricata l'icona)
iconFilename3 = 'http://maps.google.com/mapfiles/kml/paddle/ylw-blank-lv.png';

% Vettore di caratteri bianchi per non fare comparire il nome di fianco alle icone nel file kml
Name_blank = blanks(length(E(:,5))); 

filename1 = [working_path '\' fireball_name '_Trajectory_' num2str(Nstaz) '_stations.kml'];
kmlwritepoint(filename1, E(:, 5), E(:, 6), 'Icon', iconFilename3, 'IconScale', 0.5, 'Name', Name_blank);
% Alternativa: disegna una linea gialla continua al suolo
% kmlwriteline(filename1, E(:, 5), E(:, 6), 'Color','yellow','Width', 10);

filename2 = [working_path '\' fireball_name '_Trajectory_and_height_' num2str(Nstaz) '_stations.kml'];
kmlwritepoint(filename2, E(:, 5), E(:, 6), 1000*E(:, 7), 'Icon', iconFilename1, 'IconScale', 1, 'Name', Name_blank);
% Alternativa: usa i marker rossi per disegnare i punti della traiettoria
% kmlwritepoint(filename2, E(:, 5), E(:, 6), 1000*E(:, 7), 'Color','red', 'Name', Name_blank);

% Legge filename2 e aggiunge la linea che arriva al suolo nel file kml della traiettoria in
% quota
filename3 = [working_path '\' fireball_name '_Trajectory_and_height_' num2str(Nstaz) '_stations_B.kml'];
fidi=fopen(filename2,'r');
fido=fopen(filename3,'w');
while ~feof(fidi)
  l=fgetl(fidi);   % read line
  if strfind(l,'</Point>') % Stringa da sostituire
    % modify line here
    l='<extrude>1</extrude></Point>';
  end
  fprintf(fido,'%s',l);  % 'fgetl returns \n so it's embedded
end
fidi=fclose(fidi);
fido=fclose(fido);

% Open file kml della traiettoria al suolo in Google Earth per Windows
% winopen(filename1) 

% Comandi per Linux
% cmd = 'googleearth ';
% fullfilename = fullfile(pwd, filename);   
% system([cmd fullfilename])

disp('   ')
disp(' END TRAJECTORY COMPUTATION  ')
disp('   ')

end