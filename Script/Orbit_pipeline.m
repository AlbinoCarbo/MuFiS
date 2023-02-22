% ORBIT FOR FIREBALL PIPELINE
%
% Script per il calcolo dell'orbita eliocentrica di un bolide nella pipeline di PRISMA.
% Corregge per l'attrazione e la rotazione terrestre usando la classica tecnica della "attrazione zenitale".
% Usa il metodo Monte Carlo per calcolare tutte le possibili orbite compatibili con le osservazioni.
% L'algoritmo è tratto da: Z. Ceplecha, Bull. Astron. Inst. Czechosl. 38 (1987), 222-234.
%
% NOTA1: Nella generazione dei cloni Monte Carlo c'è un ciclo if che
% impedisce alla velocità all'infinito di scendere al di sotto degli 11.4
% km/s. Se così fosse l'orbita non potrebbe più essere eliocentrica e non
% sarebbe uno scenario compatibile con quanto osservato.
%
% NOTA2: Per un metodo più accurato e completo per il calcolo dell'orbita, che tenga conto dello 
% schiacciamento terrestre e della gravità di Luna, Sole e pianeti vedi:
%
% Vasily Dmitriev, Valery Lupovka, Maria Gritsevich, "Orbit determination
% based on meteor observations using numerical integration of equations of
% motion", Planetary and Space Science, Volume 117, November 2015, Pages
% 223-235.
% 
% FUNZIONI ESTERNE RICHIESTE:
%
% sunlong.m, calcola la longitudine eclittica del Sole in radianti (J2000.0)
%
% distanzaterrasole.m, calcola la distanza eliocentrica Terra-Sole in U.A.
%
% temposiderale.m, funzione per il calcolo del tempo siderale medio locale in radianti al momento della
% osservazione del fireball.
%
% EQAzimut.m, funzione che trasforma da coordinate equatoriali alla data in azimutali con azimut
% contato da nord verso est.
%
% Azimut2EQ.m, funzione per il passagio da coordinate azimutali (azimut contato da nord verso
% est) a coordinate equatoriali alla data.
%
% EQ2EQ2000.m, funzione per il passaggio dalla coordinate equatoriali
% all'equinozio medio della data all'equinozio J2000.0
%
% meccanica_celeste_eliocentrica.m, funzione per il calcolo degli elementi orbitali
% eliocentrici a partire dalle coordinate equatoriali geocentriche vere del radiante, la
% velocità geocentrica vera, la longitudine eclittica geocentrica del Sole e la distanza
% Terra-Sole.
%
% meccanica_celeste_baricentrica.m, funzione per il calcolo degli elementi orbitali
% baricentrici (ossia rispetto al baricentro del Sistema Solare), a partire dalle coordinate equatoriali 
% geocentriche vere del radiante, la velocità geocentrica vera, la longitudine eclittica geocentrica 
% del Sole e la distanza Terra-baricentro del Sistema Solare.
%
% identification.m, funzione per il calcolo di U, theta e fi per l'identificazione 
% della origine del meteoroide secondo Valsecchi et al., Mon. Not. R. Astron. Soc. 304, 743±750 (1999).
%
% Orbit_plotter_monte_carlo.m, calcola, plotta e salva una figura con le orbite - in scala e con il corretto orientamento reciproco - di tutti i pianeti e del pianeta nano Plutone. 
% Alle orbite planetarie viene sovrapposta l'orbita nominale del meteoroide progenitore del bolide più quelle (eventuali, dipende dall'input) 
% dei cloni Monte Carlo. Fa un plot a punti random anche della Fascia Principale (fra 2.2 e 3.6 UA) e della Fascia di Kuiper (fra 30 e 55 UA).
%
% outliers_deleter.m, funzione per il calcolo degli elementi di un vettore
% distanti più di 3 deviazioni standard dalla mediana del vettore stesso.
% Utile per cancellari scenari Monte Carlo troppo anomali.
%
% INPUT:
% working_path = path della cartella di salvataggio risultati
% fireball_name = nome del fireball
% AR_Rad = ascensione retta alla data del radiante apparente geocentrico del fireball (radianti)
% Dec_Rad = declinazione alla data del radiante apparente geocentrico del fireball (radianti)
% incertezza_radiante = incertezza angolare del radiante in gradi (questo valore arriva dalla triangolazione con N stazioni)
% Long_fin = Longitudine del punto osservato finale in gradi
% Lat_fin = Latitudine del punto osservato finale in gradi
% V_inf = velocità all'infinito del meteoroide in km/s
% SV_inf = incertezza della velocità all'infinito in km/s (dal modello
% dinamico del bolide)
% quota_fin = quota del punto finale in km
% sQuota_fin = incertezza quota finale del bolide in km (dal modello
% dinamico del bolide)
% N = numero di scenari Monte Carlo di cui calcolare l'orbita del
% meteoroide progenitore. In questo modo si ottiene un set di orbite tutte
% compatibili con il fireball osservato e, visualizzandolecon Celestia, si ha una
% percezione visiva delle incertezze relative all'origine del meteoroide.
% 
% OUTPUT: 
% salvataggio su file delle coordinate vere del radiante, della velocità geocentrica ed eliocentrica/baricentrica vera del meteoroide 
% al momento dell'incontro con la Terra, dell'invariante di Tisserand e degli elementi orbitali eliocentrici/baricentrici - nominali - del meteoroide. 
% Le incertezze associate agli elementi orbitali sono ottenute con il metodo Monte Carlo (default 100 cloni compatibili con le osservazioni).
%
% Salvataggio su file .ssc sia degli elementi orbitali eliocentrici nominali sia di tutti gli elementi orbitali dei meteoroidi
% progenitori virtuali compatibili con le osservazioni. Questi due file possono essere visualizzati con il planetario 3D Celestia.
% 
% Salvataggio di due figure in formato .bmp con il plot sia dell'orbita eliocentrica nominale sia di quelle con le orbite dei cloni Monte Carlo. Le figure
% riportano le posizioni dei pianeti al momento della caduta del meteoroide.
%
% Albino Carbognani (OAVdA)
%
% Versione del 9 aprile 2020

function []=Orbit_pipeline(working_path, fireball_name, jd0, UT, AR_Rad, Dec_Rad, incertezza_radiante, Lat_fin, Long_fin, V_inf, SV_inf, quota_fin, Squota_fin, N)

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('             ORBIT FOR FIREBALL (MONTE CARLO ROUTINE)                  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')

% Costanti

C=57.29577951;                     % Conversione da gradi a radianti e viceversa
epsilon=23.4392911/C;              % Obliquità dell'eclittica (radianti)
L_sun=sunlong(jd0+UT/24);          % Longitudine eclittica del Sole in radianti (J2000.0)
RE=distanzaterrasole(jd0+UT/24);   % Distanza eliocentrica Terra-Sole in UA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dati di input sul meteoroide residuo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' READ INPUT DATA WITH UNCERTAINTIES ')
disp('   ')

AR_R0=AR_Rad;                      % AR del radiante apparente in gradi ma trasformato in radianti (alla data)
sAR_R=incertezza_radiante/C;       % Incertezza AR del radiante apparente (radianti). 

Dec_R0=Dec_Rad;                    % Dec del radiante apparente in gradi ma trasformato in radianti (alla data)
sDec_R=incertezza_radiante/C;      % Incertezza Dec del radiante apparente (radianti).

V_inf0=V_inf;                      % Velocità all'infinito in km/s
sV_inf=SV_inf;                     % Incertezza velocità all'infinito in km/s (dal modello dinamico del bolide)

Long_fin0=Long_fin/C;              % Longitudine del punto finale osservato del bolide (radianti)
sLong_fin=0.001/C;                 % Incertezza longitudine punto finale (radianti). VALORE FISSO.

Lat_fin0=Lat_fin/C;                % Latitudine del punto finale osservato del bolide (radianti)
sLat_fin=0.001/C;                  % Incertezza latitudine punto finale (radianti). VALORE FISSO.

quota_fin0=quota_fin;              % Quota del punto finale del bolide in km (dal modello dinamico del bolide)
sQuota_fin=Squota_fin;             % Incertezza quota finale del bolide in km (dal modello dinamico del bolide)

N=N+1;                             % Viene aumentato di 1 il numero di cloni Monte Carlo perché il primo è il meteoroide nominale

%*********FINE DELLA FASE DI INPUT ED INIZIO DEI CALCOLI*******************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparazione per la routine Monte Carlo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocazione dei vettori di output delle coordinate del radiante corrette per velocità e rotazione terrestre 
ARr=zeros(1, N); 
DECr=zeros(1, N); 
AR_g=zeros(1, N);
DEC_g=zeros(1, N);
AR2000_g=zeros(1, N); 
DEC2000_g=zeros(1, N);
Lr=zeros(1, N); 
Br=zeros(1, N);
ARr_B=zeros(1, N); 
DECr_B=zeros(1, N);
Lr_B=zeros(1, N); 
Br_B=zeros(1, N);

% Pre-allocazione del vettore della velocità geocentrica
Vg=zeros(1, N);

% Pre-allocazione dei vettori di output sugli elementi orbitali
% eliocentrici
A=zeros(1, N);       % Semiasse maggiore
I=zeros(1, N);       % Inclinazione orbitale
NA=zeros(1, N);      % Longitudine del nodo ascendente
E=zeros(1, N);       % Eccentricità
AP=zeros(1, N);      % Longitudine del perielio
T=zeros(1, N);       % JD passaggio al perielio
V_helio=zeros(1, N); % Velocità eliocentrica
vera=zeros(1, N);    % Anomalia vera

% Pre-allocazione dei vettori di output sugli elementi orbitali
% baricentrici
A_B=zeros(1, N);    % Semiasse maggiore
I_B=zeros(1, N);    % Inclinazione orbitale
NA_B=zeros(1, N);   % Longitudine del nodo ascendente
E_B=zeros(1, N);    % Eccentricità
AP_B=zeros(1, N);   % Longitudine del perielio
T_B=zeros(1, N);    % JD passaggio al perielio
V_bari=zeros(1, N); % Velocità eliocentrica
vera_B=zeros(1, N); % Anomalia vera

% Pre-allocazione dei vettori di output per l'identificazione
% dell'eventuale corpo progenitore del meteoroide originario

U=zeros(1, N); 
Ux=zeros(1, N); 
Uy=zeros(1, N); 
Uz=zeros(1, N);
theta=zeros(1, N); 
fi=zeros(1, N);
AA=zeros(1, N);
EE=zeros(1, N);
II=zeros(1, N);

disp('   ')
disp(' START MONTE CARLO COMPUTATION FOR HELIOCENTRIC/BARYCENTRIC ORBIT  ')
disp('   ')

for i=1:N
    
% NOTA: con l'istruzione if che segue la prima orbita calcolata è quella nominale.
% Il computo Monte Carlo serve per la stima della incertezza da assegnare
% ai vari elementi orbitali nominali.
    
    if i==1
        orbita_nominale=0;
    else
        orbita_nominale=1;
    end
    
% Generazione dello scenario Monte Carlo compatibile con le osservazioni
% e con un'orbita eliocentrica.

AR_R = AR_R0 + orbita_nominale*(sAR_R)*(randn);                 % AR del radiante apparente (radianti)
Dec_R = Dec_R0 + orbita_nominale*(sDec_R)*(randn);              % Dec del radiante apparente (radianti)
V_inf = V_inf0 + orbita_nominale*(sV_inf)*(randn);              % Velocità all'infinito (km/s)

if V_inf < 11.4
V_inf = V_inf0 + orbita_nominale*(sV_inf)*(randn);              % Estrazione nuova velocità all'infinito se V_inf è inferiore agli 11.4 km/s (km/s)    
end                                                             % Nota: se V_inf < 11.4 km/s il meteoroide non può avere un'orbita eliocentrica!

Long_fin = Long_fin0 + orbita_nominale*(sLong_fin)*(randn);     % Longitudine del punto finale osservato del bolide (radianti)
Lat_fin = Lat_fin0 + orbita_nominale*(sLat_fin)*(randn);        % Latitudine del punto finale osservato del bolide (radianti)
quota_fin = quota_fin0 + orbita_nominale*(sQuota_fin)*(randn);  % Quota del punto finale del bolide (km)

% NOTA: randn è un numero casuale con distribuzione normale standard, ossia con valore medio
% zero e deviazione standard 1. Lo scenario Monte Carlo compatibile con le
% osservazioni viene generato per i > 2. Per i = 1 si ha l'orbita nominale.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocità di rotazione della Terra alla latitudine del punto finale osservato (km/s) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=geocradius(C*Lat_fin, 'WGS84')/1000;           % Raggio geocentrico del punto finale del bolide (km)

VE=2*pi*(R+quota_fin)*cos(Lat_fin)/86164.09;     % Velocità di rotazione della Terra nel punto finale del bolide (km/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo della velocità geocentrica all'infinito %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Versori geocentrici del radiante apparente medio
csi=cos(Dec_R).*cos(AR_R);
eta=cos(Dec_R).*sin(AR_R);
zeta=sin(Dec_R);

% Componenti geocentriche della velocità all'infinito in km/s
Vx=V_inf*csi; Vy=V_inf*eta; Vz=V_inf*zeta;  

% Componenti della velocità geocentrica all'infinito corrette per la
% rotazione terrestre

[theta_E] = temposiderale(jd0, UT, Long_fin); % Tempo siderale nel punto finale osservato (in radianti)

Vxc=Vx-VE*cos(theta_E+pi/2);   % Attenzione: l'argomento del coseno è il tempo siderale del punto est
Vyc=Vy-VE*sin(theta_E+pi/2);   % Attenzione: l'argomento del seno è il tempo siderale del punto est
Vzc=Vz;

% Modulo della velocità geocentrica all'infinito corretta per la rotazione
% terrestre
Vc=sqrt(Vxc^2+Vyc^2+Vzc^2);

% Modulo della velocità geocentrica all'infinito corretta anche per l'attrazione
% terrestre (km/s)

Vg(i)=sqrt(Vc^2-797201.0/(R+quota_fin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correzione del radiante apparente per la gravità terrestre %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ascensione retta del radiante corretto per la rotazione terrestre

X=Vxc/Vc; Y=Vyc/Vc;

AR_c=atan2(Y, X);     % Ascensione retta geocentrica del radiante corretto per la rotazione terrestre (compresa fra -pi e +pi)

if AR_c < 0
    AR_c=AR_c+2*pi;   % Riduzione dell'AR ad un valore sempre positivo compreso fra 0 e 2*pi
end

% Declinazione del radiante corretto per la rotazione terrestre

Delta_c=asin(Vzc/Vc);

% Passaggio a coordinate azimutali

[Hc, Ac] = EQAzimut(Lat_fin, AR_c, Delta_c, theta_E);

% Calcolo della distanza zenitale:

Zc=pi/2-Hc;

Delta_Zc=2*atan((Vc-Vg(i))*tan(Zc/2)/(Vc+Vg(i)));

Zg=Zc+Delta_Zc; % Distanza zenitale del radiante corretto per la gravità terrestre

% Calcolo dell'azimut del radiante corretto per la gravità terrestre

Ag=Ac;

% Coordinate equatoriali del radiante corretto per la gravità terrestre (in
% radianti, equinozio alla data)

[AR_g(i), DEC_g(i)] = Azimut2EQ(Lat_fin, Ag, pi/2-Zg, theta_E);

% Coordinate equatoriali del radiante corretto per la gravità terrestre (in
% radianti, equinozio J2000.0)

[AR2000_g(i), DEC2000_g(i)] = EQ2EQ2000(AR_g(i), DEC_g(i), jd0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcolo dell'orbita eliocentrica del meteoroide     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Lr(i), Br(i), ARr(i), DECr(i), A(i), I(i), NA(i), E(i), AP(i), T(i), V_helio(i), vera(i)]=meccanica_celeste_eliocentrica(jd0+UT/24, AR2000_g(i), DEC2000_g(i), Vg(i), L_sun, RE);

% NOTA:
% L_sun è la longitudine eclittica del Sole in radianti (J2000.0)
% RE è la distanza eliocentrica Terra-Sole in UA

% Trasformazione della velocità eliocentrica da UA/d a km/s

V_helio(i)=(1731.456829)*V_helio(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcolo dell'orbita baricentrica del meteoroide     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Lr_B(i), Br_B(i), ARr_B(i), DECr_B(i), A_B(i), I_B(i), NA_B(i), E_B(i), AP_B(i), T_B(i), V_bari(i), vera_B(i)]=meccanica_celeste_baricentrica(jd0+UT/24, AR2000_g(i), DEC2000_g(i), Vg(i), L_sun, RE);

% Trasformazione della velocità baricentrica da UA/d a km/s

V_bari(i)=(1731.456829)*V_bari(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Identificazione del corpo progenitore del meteoroide       %
%(Valsecchi et al., Mon. Not. R. Astron. Soc. 304, 743-750 (1999)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[U(i), Ux(i), Uy(i), Uz(i), theta(i), fi(i), AA(i), EE(i), II(i)]=identification(AR2000_g(i), DEC2000_g(i), Vg(i), epsilon, L_sun-pi);
% NOTA:
% epsilon = obliquità dell'Eclittica
% L_sun-pi = longitudine Eclittica della Terra


end 

disp('   ')
disp(' END MONTE CARLO COMPUTATION  ')
disp('   ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del vettore dell'invariante di Tisserand rispetto a Giove %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE TISSERAND INVARIANT  ')
disp('   ')

ap=5.20336301; % Semiasse maggiore dell'orbita di Giove in UA

TJ=(ap./A)+2*(sqrt((A/ap).*(1-E.*E))).*cos(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute probability for heliocentric hyperbolic orbit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyp_counter=0;

for i=1:N
    
    if E(i)>1
        
        hyp_counter=hyp_counter+1;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute probability for baricentric hyperbolic orbit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyp_counter_B=0;

for i=1:N
    
    if E_B(i)>1
        
        hyp_counter_B=hyp_counter_B+1;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Salvataggio dei risultati in un file di testo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' SAVE HELIOCENTRIC/BARYCENTRIC NOMINAL ORBIT DATA IN A .TXT FILE  ')
disp('   ')
fid = fopen([working_path '\' fireball_name '_Orbit.txt'],'w'); % Path di salvataggio del file di dati
fprintf(fid, ['Fireball: ' fireball_name '\n']);
fprintf(fid, 'JD0: %10.6f \n', jd0);
fprintf(fid, 'UT:  %10.6f \n', UT);
fprintf(fid, '    \n');
fprintf(fid, '**************************************************************************************************************************************************************   \n');
fprintf(fid, 'Starting values and uncertainties to which the Monte Carlo method is applied to propagate them to the true radiant and to the heliocentric orbital elements\n');
fprintf(fid, '    \n');
fprintf(fid, 'RA of the apparent geocentric radiant (degrees)              %10.4f \n', C*AR_R0);
fprintf(fid, 'RA uncertainty of the apparent geocentric radiant (degrees)  %10.4f \n', incertezza_radiante);
fprintf(fid, 'Dec of the apparent geocentric radiant (degrees)             %10.4f \n', C*Dec_R0);
fprintf(fid, 'Dec uncertainty of the apparent geocentric radiant (degrees) %10.4f \n', incertezza_radiante);
fprintf(fid, 'Pre-atmospheric velocity (km/s)                              %10.4f \n', V_inf0);
fprintf(fid, 'Pre-atmospheric velocity uncertainty (km/s)                  %10.4f \n', sV_inf);
fprintf(fid, 'Final point longitude (degree)                               %10.4f \n', C*Long_fin0);
fprintf(fid, 'Final point longitude uncertainty (degree)                   %10.4f \n', sLong_fin*C);
fprintf(fid, 'Final point latitude (degree)                                %10.4f \n', C*Lat_fin0);
fprintf(fid, 'Final point latitude uncertainty (degree)                    %10.4f \n', sLat_fin*C);
fprintf(fid, 'Final height (km)                                            %10.4f \n', quota_fin0);
fprintf(fid, 'Final height uncertainty (km)                                %10.4f \n', sQuota_fin);
fprintf(fid, '***************************************************************************************************************************************************************   \n');
fprintf(fid, '    \n');
fprintf(fid, 'CELESTIAL COORDINATES OF THE TRUE RADIANT \n');
fprintf(fid, '    \n');
fprintf(fid, 'Geocentric equatorial coordinates of the true radian (mean equinox to date) \n');
fprintf(fid, '  AR (°) \t\t\t\t\t Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', AR_g(1)*C, std(outliers_deleter(AR_g*C)), DEC_g(1)*C, std(outliers_deleter(DEC_g*C)));
fprintf(fid, '    \n');
fprintf(fid, 'Geocentric equatorial coordinates of the true radiant (equinox J2000.0) \n');
fprintf(fid, '  RA (°) \t\t\t\t\t Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', AR2000_g(1)*C, std(outliers_deleter(AR2000_g*C)), DEC2000_g(1)*C, std(outliers_deleter(DEC2000_g*C)));
fprintf(fid, '    \n');
fprintf(fid, 'Heliocentric equatorial coordinates of the true radiant (equinox J2000.0) \n');
fprintf(fid, '  RA (°)  \t\t\t\t\t   Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', ARr(1)*C, std(outliers_deleter(ARr*C)), DECr(1)*C, std(outliers_deleter(DECr*C)));
fprintf(fid, '    \n');
fprintf(fid, 'Heliocentric ecliptic coordinates of the true radiant (equinox J2000.0) \n');
fprintf(fid, '  RA (°)  \t\t\t\t\t   Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', Lr(1)*C, std(outliers_deleter(Lr*C)), Br(1)*C, std(outliers_deleter(Br*C)));
fprintf(fid, '    \n');
fprintf(fid, 'Barycentric equatorial coordinates of the true radiant (equinox J2000.0) \n');
fprintf(fid, '  RA (°)  \t\t\t\t\t   Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', ARr_B(1)*C, std(outliers_deleter(ARr_B*C)), DECr_B(1)*C, std(outliers_deleter(DECr_B*C)));
fprintf(fid, '    \n');
fprintf(fid, 'Barycentric ecliptic coordinates of the true radiant (equinox J2000.0) \n');
fprintf(fid, '  RA (°)  \t\t\t\t\t   Dec (°)\n');
fprintf(fid, '    \n');
fprintf(fid, '%10.2f \t pm %10.2f %10.2f \t pm %10.2f \n', Lr_B(1)*C, std(outliers_deleter(Lr_B*C)), Br_B(1)*C, std(outliers_deleter(Br_B*C)));

fprintf(fid, '    \n');
fprintf(fid, '    \n');
fprintf(fid, 'VALUES OF THE TRUE GEOCENTRIC VELOCITY, TRUE HELIOCENTRIC/BARYCENTRIC VELOCITY AND JUPITER TISSERAND INVARIANT\n');
fprintf(fid, '    \n');
fprintf(fid, 'True geocentric speed (km/s): %10.2f \t pm %10.2f \n', Vg(1), std(Vg));
fprintf(fid, '    \n');
fprintf(fid, 'True heliocentric speed (km/s): %10.2f \t pm %10.2f \n', V_helio(1), std(V_helio));
fprintf(fid, '    \n');
fprintf(fid, 'True barycentric speed (km/s): %10.2f \t pm %10.2f \n', V_bari(1), std(V_bari));
fprintf(fid, '    \n');
fprintf(fid, 'Jupiter Tisserand Invariant: %10.2f \t pm %10.2f \n', TJ(1), std(TJ));
fprintf(fid, '    \n');
fprintf(fid, 'Probability for heliocentric hyperbolic orbit (per cent):    %10.2f \n', 100*hyp_counter/N);
fprintf(fid, '    \n');
fprintf(fid, 'Probability for barycentric hyperbolic orbit (per cent):     %10.2f \n', 100*hyp_counter_B/N);
fprintf(fid, '    \n');
fprintf(fid, 'METEOROID NOMINAL HELIOCENTRIC ORBITAL ELEMENTS (J2000.0)\n');
fprintf(fid, 'Uncertainties are obtained with the Monte Carlo method\n');
fprintf(fid, 'Clones number = %10.0f    \n', N);
fprintf(fid, '    \n');
fprintf(fid, 'Semi-major axis (AU)               %10.3f \t pm %10.3f \n', A(1), std(A));

if E(1) < 1
fprintf(fid, 'Orbital period (years)             %10.3f \t pm %10.3f \n', sqrt(A(1)^3), std(sqrt(A.^3)));
else
fprintf(fid, 'Orbital period (years)             infty \n');   % Per orbite paraboliche o iperboliche 
end

% Riduzione longitudine del perielio eliocentrico fra 0 e 2*pi
LOPh=AP(1)+NA(1);
if LOPh>2*pi
    LOPh=LOPh-2*pi;
end

fprintf(fid, 'Orbit inclination (°)              %10.3f \t pm %10.3f \n', I(1)*C, std(I*C));
fprintf(fid, 'Longitude ascending node (°)       %10.5f \t pm %10.3f \n', NA(1)*C, std(NA*C));
fprintf(fid, 'Eccentricity                       %10.3f \t pm %10.3f \n', E(1), std(E));
fprintf(fid, 'Argument of perihelion (°)         %10.3f \t pm %10.3f \n', AP(1)*C, std(AP*C));
fprintf(fid, 'Longitude of perihelion (°)        %10.3f \t pm %10.3f \n', (LOPh)*C, sqrt((std(AP*C))^2+(std(NA*C))^2));
fprintf(fid, 'True anomaly (°)                   %10.3f \t pm %10.3f \n', vera(1)*C, std(vera*C));
fprintf(fid, 'Perihelion passage (JD)            %10.2f \t pm %10.0f \n', T(1), std(T));
fprintf(fid, 'Perihelion distance (AU)           %10.3f \t pm %10.3f \n', A(1)*(1-E(1)), std(A.*(1-E)));

if E(1) < 1
fprintf(fid, 'Aphelion distance (AU)             %10.3f \t pm %10.3f \n', A(1)*(1+E(1)), std(A.*(1+E)));
else
fprintf(fid, 'Aphelion distance (AU)             infty \n');  % Per orbite paraboliche o iperboliche   
end

fprintf(fid, '    \n');
fprintf(fid, '    \n');
fprintf(fid, 'METEOROID NOMINAL BARYCENTRIC ORBITAL ELEMENTS (J2000.0)\n');
fprintf(fid, 'Uncertainties are obtained with the Monte Carlo method\n');
fprintf(fid, '    \n');
fprintf(fid, 'Semi-major axis (AU)               %10.3f \t pm %10.3f \n', A_B(1), std(A_B));

if E_B(1) < 1
fprintf(fid, 'Orbital period (years)             %10.3f \t pm %10.3f \n', sqrt(A_B(1)^3), std(sqrt(A_B.^3)));
else
fprintf(fid, 'Orbital period (years)             infty \n');   % Per orbite paraboliche o iperboliche 
end

% Riduzione longitudine del perielio baricentrico fra 0 e 2*pi
LOPb=AP_B(1)+NA_B(1);
if LOPb > 2*pi
    LOPb=LOPb-2*pi;
end

fprintf(fid, 'Orbit inclination (°)              %10.3f \t pm %10.3f \n', I_B(1)*C, std(I_B*C));
fprintf(fid, 'Longitude ascending node (°)       %10.5f \t pm %10.3f \n', NA_B(1)*C, std(NA_B*C));
fprintf(fid, 'Eccentricity                       %10.3f \t pm %10.3f \n', E_B(1), std(E_B));
fprintf(fid, 'Argument of perihelion (°)         %10.3f \t pm %10.3f \n', AP_B(1)*C, std(AP_B*C));
fprintf(fid, 'Longitude of perihelion (°)        %10.3f \t pm %10.3f \n', (LOPb)*C, sqrt((std(AP_B*C))^2+(std(NA_B*C))^2));
fprintf(fid, 'True anomaly (°)                   %10.3f \t pm %10.3f \n', vera_B(1)*C, std(vera_B*C));
fprintf(fid, 'Perihelion passage (JD)            %10.2f \t pm %10.0f \n', T_B(1), std(T_B));
fprintf(fid, 'Perihelion distance (AU)           %10.3f \t pm %10.3f \n', A_B(1)*(1-E_B(1)), std(A_B.*(1-E_B)));

if E_B(1) < 1
fprintf(fid, 'Aphelion distance (AU)             %10.3f \t pm %10.3f \n', A_B(1)*(1+E_B(1)), std(A_B.*(1+E_B)));
else
fprintf(fid, 'Aphelion distance (AU)             infty \n');  % Per orbite paraboliche o iperboliche   
end

fprintf(fid, '    \n');
fprintf(fid, '    \n');
fprintf(fid, 'IDENTIFICATION DATA FOR METEOROID PROGENITOR (Valsecchi et al., Mon. Not. R. Astron. Soc. 304, 743-750 (1999))\n');
fprintf(fid, 'Uncertainties are obtained with the Monte Carlo method\n');
fprintf(fid, 'Clones number = %10.0f    \n', N);
fprintf(fid, '    \n');
fprintf(fid, 'Geocentric velocity Ux (in units of the terrestrial one) %10.3f \t pm %10.2f \n', Ux(1), std(Ux));
fprintf(fid, 'Geocentric velocity Uy (in units of the terrestrial one) %10.3f \t pm %10.2f \n', Uy(1), std(Uy));
fprintf(fid, 'Geocentric velocity Uz (in units of the terrestrial one) %10.3f \t pm %10.2f \n', Uz(1), std(Uz));
fprintf(fid, '    \n');
fprintf(fid, 'Geocentric velocity U (in units of the terrestrial one)  %10.3f \t pm %10.2f \n', U(1), std(U));
fprintf(fid, 'Theta (°)                                                %10.3f \t pm %10.2f \n', theta(1)*C, std(theta*C));
fprintf(fid, 'Fi (°)                                                   %10.3f \t pm %10.2f \n', fi(1)*C, std(fi*C));
fprintf(fid, '    \n');
fprintf(fid, 'Orbital elements from geocentric velocity U\n');
fprintf(fid, 'Semi-major axis (AU)                                     %10.3f \t pm %10.3f \n', AA(1), std(AA));
fprintf(fid, 'Eccentricity                                             %10.3f \t pm %10.3f \n', EE(1), std(EE));
fprintf(fid, 'Orbit inclination (°)                                    %10.3f \t pm %10.3f \n', II(1)*C, std(II));

% Chiusura del file di output
status = fclose(fid);
      
disp('   ')
disp(' PLOT PROGENITOR METEOROID NOMINAL HELIOCENTRIC ORBIT IN SOLAR SYSTEM ')
disp('   ')

% Plot dell'orbita eliocentrica nominale del fireball
Orbit_plotter_monte_carlo(fireball_name, working_path, jd0+UT/24, A(1), E(1), (AP(1)+NA(1))*C, 1);

disp('   ')
disp(' SAVE NOMINAL HELIOCENTRIC ORBIT IN A .SSC FILE VIEWABLE WITH CELESTIA ')
disp('   ')
fid0 = fopen([working_path '\' fireball_name '_Nominal_Orbit_Celestia.ssc'],'w'); % Path di salvataggio del file di dati
fprintf(fid0, ['"Fireball ' fireball_name ' " "Sol" \n']);
fprintf(fid0, ' { \n');
fprintf(fid0, '	Class "planet"  \n');            % La classe non è asteroid ma planet così i colori delle orbite visualizzate da Celestia sono diversi
fprintf(fid0, '	Texture "asteroid.jpg"     \n'); 
fprintf(fid0, '	Color [1.000 0.945 0.881]  \n'); % Colore dell'asteroide
fprintf(fid0, '	BlendTexture true          \n');
fprintf(fid0, '	Mesh "asteroid.cms"        \n');
fprintf(fid0, '	Radius          38.94      \n');
fprintf(fid0, '      \n');
fprintf(fid0, '	EllipticalOrbit \n');
fprintf(fid0, '	{ \n');
fprintf(fid0, '	Epoch           %10.4f \n', T(1));
fprintf(fid0, '	Period          %10.4f \n', sqrt(A(1)^3));
fprintf(fid0, '	SemiMajorAxis   %10.4f \n', A(1)); 
fprintf(fid0, '	Eccentricity    %10.4f \n', E(1)); 
fprintf(fid0, '	Inclination     %10.4f \n', I(1)*C);
fprintf(fid0, '	AscendingNode   %10.4f \n', NA(1)*C);
fprintf(fid0, '	ArgOfPericenter %10.4f \n', AP(1)*C);
fprintf(fid0, '	MeanAnomaly     %10.4f \n', 0.0);
fprintf(fid0, '	} \n');
fprintf(fid0, ' } \n');

status = fclose(fid0);

disp('   ')
disp(' SAVE MONTE CARLO CLONE HELIOCENTRIC ORBITS IN A .SSC FILE VIEWABLE WITH CELESTIA ')
disp('   ')

% Salvataggio delle possibili orbite Monte Carlo in un file .ssc visualizzabile con Celestia
fid1 = fopen([working_path '\' fireball_name '_Orbits_Monte_Carlo_Celestia.ssc'],'w'); % Path di salvataggio del file di dati

for i=2:N     % Si inizia con i = 2 perché per i = 1 si ha l'orbita eliocentrica nominale del meteoroide progenitore salvata nel file .ssc precedente
     
numero = num2str(i); % Conversione da numero a stringa
fprintf(fid1, ['"Fireball MC ' fireball_name '  ' numero ' " "Sol" \n']);
fprintf(fid1, ' { \n');
fprintf(fid1, '	Class "asteroid"  \n');
fprintf(fid1, '	Texture "asteroid.jpg"     \n'); 
fprintf(fid1, '	Color [1.000 0.945 0.881]  \n');
fprintf(fid1, '	BlendTexture true          \n');
fprintf(fid1, '	Mesh "asteroid.cms"        \n');
fprintf(fid1, '	Radius          38.94      \n');
fprintf(fid1, '      \n');
fprintf(fid1, '	EllipticalOrbit \n');
fprintf(fid1, '	{ \n');
fprintf(fid1, '	Epoch           %10.4f \n', T(i));
fprintf(fid1, '	Period          %10.4f \n', sqrt(A(i)^3));
fprintf(fid1, '	SemiMajorAxis   %10.4f \n', A(i)); 
fprintf(fid1, '	Eccentricity    %10.4f \n', E(i)); 
fprintf(fid1, '	Inclination     %10.4f \n', I(i)*C);
fprintf(fid1, '	AscendingNode   %10.4f \n', NA(i)*C);
fprintf(fid1, '	ArgOfPericenter %10.4f \n', AP(i)*C);
fprintf(fid1, '	MeanAnomaly     %10.4f \n', 0.0);
fprintf(fid1, '	} \n');
fprintf(fid1, ' } \n');
fprintf(fid1, '   \n');

end

% Chiusura del file di output .ssc per i cloni
status = fclose(fid1);
          
disp('   ') 
disp(' PLOT MONTE CARLO PROGENITOR METEOROIDS HELIOCENTRIC ORBITS IN SOLAR SYSTEM ')
disp('   ')

% Plot delle orbite eliocentriche di tutti i cloni Monte Carlo e di quella nominale
Orbit_plotter_monte_carlo(fireball_name, working_path, jd0+UT/24, A, E, (AP+NA)*C, N);

disp('   ')
disp(' END HELIOCENTRIC/BARYCENTRIC ORBITS COMPUTATION  ')
disp('   ')

end

