% METEOROID DYNAMIC
%
% Script di Matlab per la risoluzione delle equazioni del moto di un fireball e 
% la determinazione dei principali parametri del meteoroide: coefficiente di drag (fissato), coefficiente di ablazione
% rapporto M/S all'infinito e velocità all'infinito. 
%
% Riferimento per le equazioni del moto: Kalenichenko, A&A 448,(2007).
% Lo script fitta anche il modello di Gritsevich (Solar System Research, 2007) per ottenere i parametri
% alpha e beta, utili per capire se c'è stata una caduta di meteoriti o meno. Viene anche ricostruita, usando 
% i risultati del modello dinamico, la curva di luce teorica del fireball, ossia magnitudine assoluta in funzione del
% tempo. L'efficenza luminosa tau è fissata a 0.04.
%
% Lo script calcola anche quota, velocità e accelerazione per il punto finale della 
% traiettoria osservata usando il modello dinamico del meteoroide calcolato nel tempo finale realmente osservato. 
% Questi ultimi dati sono essenziali per il modello della fase di volo buio che fornisce la zona di ricerca delle eventuali meteoriti. 
% ATTENZIONE: questo script va usato dopo la triangolazione della traiettoria! 
%
% IPOTESI FISICHE DI BASE: 1) Il meteoroide deve mantenere la sua forma durante
% l'ablazione e non deve ruotare 2) La superfice terrestre è considerata piatta. 
% Modello di Kalenichenko, A&A 448,(2007). 
% 
% La densità dell'aria è espressa usando la legge esponenziale in modo che la densità 
% in funzione della quota si adatti al meglio possibile al modello 
% COESA del 1976 (modello non isotermico).
%
% FUNZIONI USATE: 
% binning_vector.m, funzione per il bin, con passo di 0.1 s, dei vettori
% osservati della posizione geocentrica, quota, velocità e tempi. Per
% mediare i dati nel caso di elevata rumorosità su tutta la curva.
%
% Ceplecha_model_velocity2.m, Funzione per la stima preliminare della velocità all'infinito 
% del meteoroide con il modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h).  
% Se il modello di Ceplecha non converge si passa ad un modello lineare.
%
% paramfun.m, funzione con le equazioni differenziali del moto del meteoroide. Contiene i parametri
% del fireball, ossia coefficiente di drag, coefficiente di ablazione, velocità e rapporto M/S all'infinito. 
%
% rho_atmosphere_COESA.m, Calcola la densità dell'atmosfera in funzione della quota 
% secondo il "Committee on Extension to the Standard Atmosphere" (COESA) model del 1976. 
% A questo scopo al suo interno usa la funzione "atmoscoesa" di Matlab. 
%
% best_fit_COESA_model.m, funzione per il calcolo dell'altezza di scala efficace di un modello esponenziale 
% dell'atmosfera in modo da fittare al meglio il modello COESA del 1976.
%
% monte_carlo_best_fit_parameter.m, funzione per la stima dell'incertezza 
% dei parametri di best fit del modello dinamico con il metodo Monte Carlo. 
% Genera N scenari diversi di posizione e velocità aggiungendo ai valori osservati 
% gli scarti medi dati dalla media della differenza fra best fit e osservati, 
% secondo una distribuzione gaussiana.
%
% INPUT:
% data_path = path dove si trovano i risultati della triangolazione;
% working_path = path di salvataggio dei risultati del modello dinamico;
% fireball_name = nome del fireball; 
% Nstaz = numero delle stazioni osservative; 
% Tin/Tfin, = Starting/ending time in s for fireball dynamic analysis;
% Azimut_fin = azimut della traiettoria; 
% mode = modalità di trattamento dei dati per il calcolo del modello dinamico. Sono possibili due valori:
%
%   1-mode='delete', se si vogliono cancellare dei punti all'inizio o alla fine della traiettoria osservata 
%     così si evita di fittare un modello su punti mal osservati che, tipicamente, sono quelli iniziali e finali.
%     Questa modalità è quella che fornisce i migliori risultati.
%   2-mode='smooth', se si vuole ottenere una traiettoria media binnando i
%     dati ad intervalli di 0.1 s. In questa modalità il software calcola
%     anche l'accelerazione osservata in funzione del tempo e mostra/salva le relative figure.
%
% G = valore fisso del coefficiente di drag per la fase di fireball. In
% questo modo il fit del modello fisico è facilitato perché ci sono solo 3
% costandi da determinare: ablazione, velocità e rapporto M/S all'infinito.
%
% bin_time = intervallo temporale per il bin dei dati (s).
% Nmc = Numero di scenari Monte Carlo della traiettoria osservata per il computo 
% della incertezza dei parametri di best fit del modello dinamico.
% L'incertezza fornita per ciascun parametro è la deviazione standard della media.
%
% OUTPUT: 
% Vinf = velocità all'infinito del modello dinamico (km/s) 
% sigmaVinf = incertezza velocità all'infinito del modello dinamico (km/s)
% Hfin = quota del punto finale del bolide (km) 
% sigmaHfin = incertezza della quota del punto finale del bolide (km) 
% V_minima = velocità nel punto finale del bolide (km/s) 
% final_acceleration_estimate = accelerazione nel punto finale (km/s^2)
%
% Salvataggio su file dei parametri del meteoroide che derivano dal modello dinamico e delle figure che mostrano il fit
% fra il modello del meteoroide con le velocità e le quote realmente osservate.
% Nel file vengono anche salvati il logaritmo naturale del coefficiente balistico, ln(alpha), e quello del coefficiente di
% perdita di massa, ln(beta). Questi numeri sono ottenuti dal modello dinamico e i risultati plottati con il modello di Gritsevich. 
% Questi due numeri servono per capire che probabilità ci possa essere per il ritrovamento di una meteorite (Gritsevich, ASR, 2009). 
% Nel file di testo viene salvata anche la stima teorica della magnitudine assoluta al massimo di luminosità del bolide.
% 
% ATTENZIONE1: La variabile di output "final_acceleration_estimate", che serve per la stima finale e iniziale del rapporto M/S, è l'accelerazione media della traiettoria 
% stimata con il modello empirico di Ceplecha. Se "final_acceleration_estimate" >= 0 lo script termina l'esecuzione perché si tratta di un 
% valore fisicamente impossibile che farebbe andare in tilt tutto il fit del modello dinamico del meteoroide.
%
% ATTENZIONE2: prima di iniziare l'integrazione numerica delle equazioni differenziali del moto usando la funzione "paramfun.m" viene fatta la
% verifica se esistono tempi identici misurati da due stazioni diverse. Visto che i tempi devono essere strettamente crescenti per poter 
% fare l'integrazione numerica, se ciò non accade, ossia si ha un Dt=0, viene sommato 1/5000 s in modo da superare il problema senza alterare in modo sensibile i dati osservati.
% In questo modo non si esce più dalla pipeline per evitare l'errore che restituisce "ode45" contenuto in "paramfun.m".
%
% ATTENZIONE3: alle righe 366-370 (circa) ci sono i constrain sui parametri del
% meteoroide, ossia i limiti inferiori e superiori di variazione consentita a parire dai parametri iniziali stimati. Notare che i parametri di variazione del coefficiente di drag G sono zero
% (ossia il drag è fissato) perché nelle eq. differenziali del moto c'è il rapporto G/D (D = M/S) e non si possono stimare queste due quantità separatamente.
% Quindi ho preferito fittare M/S e fissare G il cui valore è meglio stimabile.
%
% Albino Carbognani (OAS)
% Versione del 6 maggio 2020

function [Vinf, sigmaVinf, Hfin, sigmaHfin, V_minima, final_acceleration_estimate] = Meteoroid_dynamic(data_path, working_path, fireball_name, Nstaz, Tin, Tfin, Azimut_fin, mode, G, bin_time, Nmc)

% clear all

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                METEOROID DYNAMIC - Fit for single body model        %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')

% Costanti globali
C=57.29577951;                      % Conversione da gradi a radianti e viceversa
Monte_Carlo_Number=Nmc;             % Numero di cicli Monte Carlo per il computo della incertezza dei parametri di best fit del modello dinamico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lettura da file dei risultati sulla triangolazione  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_path)                                          

disp('   ')
disp(' READ TRIANGULATION RESULTS  ')
disp('   ')

s0=[fireball_name '_Trajectory_' num2str(Nstaz) '_stations_ordered.txt'];
Dat(:,:)=load(s0, '-ascii');

Tempo=Dat(:, 1);                                         % Tempo osservato dell'evento (s)
X_N=Dat(:, 2);                                           % X geocentrica in km
Y_N=Dat(:, 3);                                           % Y geocentrica in km
Z_N=Dat(:, 4);                                           % Z geocentrica in km
Lat=Dat(:, 5);                                           % Latitudine dei punti osservati della traiettoria (gradi)
Long=Dat(:, 6);                                          % Longitudine dei punti osservati della traiettoria (gradi)
Height=Dat(:, 7);                                        % Altezza del fireball in km
Velocity=Dat(:, 8);                                      % Velocità in km/s

Tempo=Tempo-Tempo(1);                                    % Conteggio del tempo a partire dal primo punto osservato che sarà Tempo(1)=0.

% INIZIO MODALITA' DELETE

tf = strcmp(mode, 'delete'); % Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
if tf == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eliminazione dei punti troppo scatterati a inizio e fine traiettoria del fireball  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' DELETE SCATTERED START/END POINTS  ')
disp('   ')

j=1; % Indice dei nuovi vettori di tempo, altezza e velocità osservati
for i=1:length(Tempo)
   
    if Tempo(i)>=Tin && Tempo(i)<=Tfin
       Tempo1(j)=Tempo(i); Height1(j)=Height(i); Velocity1(j)=Velocity(i);
       X_N1(j)=X_N(i); Y_N1(j)=Y_N(i); Z_N1(j)=Z_N(i);
       j=j+1;
    end
    
end
   
% FINE MODALITA' DELETE

else
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo della traiettoria media del fireball  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% INIZIO MODALITA' SMOOTH

disp('   ')
disp([' SMOOTH ALL TRAJECTORY POINTS IN ' num2str(bin_time) ' s BIN '])
disp('   ')

[X_N1, Y_N1, Z_N1, Height1, Velocity1, Tempo1]=binning_vector(X_N, Y_N, Z_N, Height, Velocity, Tempo, bin_time);

end

% FINE MODALITA' SMOOTH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del coseno medio della distanza zenitale della traiettoria selezionata %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocazione dei vettori per aumentare la velocità di esecuzione
L=zeros(1, length(Height1)-1);    
H=zeros(1, length(Height1)-1); 

for i=1:length(Height1)-1
L(i)=sqrt((X_N1(1)-X_N1(i+1))^2+(Y_N1(1)-Y_N1(i+1))^2+(Z_N1(1)-Z_N1(i+1))^2);
H(i)=Height1(1)-Height1(i+1);
end

COS_Zr=mean(H./L); % Coseno medio della distanza zenitale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stima preliminare della velocità all'infinito e dell'accelerazione    %
%      del meteoroide usando il modello empirico di Ceplecha (1961).      %
%  Si passa ad un modello lineare se il primo ha dei problemi di best fit.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' OBSERVED VELOCITY - CEPLECHA OR LINEAR FIT VS HEIGHT  ')
disp('   ')

Inclinazione=pi/2-acos(COS_Zr); % Inclinazione media traiettoria (radianti)

[V_infty, DV_infty, fit_velocity, Afin, A2, model_type]=Ceplecha_model_velocity2(Height, Height1, Velocity, Velocity1, Tempo1, Inclinazione);

% Accelerazione nel punto finale del bolide (km/s^2)

final_acceleration_estimate=Afin; 

if final_acceleration_estimate >= 0
    disp('   ')
    disp(' WARNING FROM METEOROID DYNAMIC FUNCTION: RAW ESTIMATE OF THE FINAL ACCELERATION ')
    disp(' IS ZERO OR POSITIVE, EXIT METEOROID DYNAMIC!')
    disp('   ')
        
    % ATTENZIONE: valori di output fittizi per la funzione "Meteoroid dynamic" in modo da proseguire
    % l'esecuzione nella "Pipeline_fireball", altrimenti da errore.
    
    Vinf=9999; 
    sigmaVinf=9999;
    Hfin=9999; 
    V_minima=9999;
    sigmaHfin=9999;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo altezza di scala e densità al suolo efficaci dell'atmosfera    %
% esponenziale in modo da fittare al meglio il modello COESA del 1976    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE BEST EXPONENTIAL ATMOSPHERIC MODEL FITTING COESA MODEL (1976)  ')
disp('   ')

[rho0_best, Hs_best] = best_fit_COESA_model(Height1);

Model_rho_atmosphere_COESA=@(x)(rho0_best*exp(1).^(-x/Hs_best));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condizioni iniziali del fireball %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=COS_Zr;                            % Coseno medio della distanza zenitale della traiettoria
V0=fit_velocity(1);                  % Velocità iniziale del fireball secondo il fit con il modello di Ceplecha o quello Lineare (km/s) <************
H0=Height1(1);                       % Quota iniziale osservata del fireball in km
rho0=Model_rho_atmosphere_COESA(H0); % Densità dell'aria alla quota iniziale osservata (kg/km^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stima dei parametri iniziali per il modello dinamico del meteoroide    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantità ausiliarie
Hf=Height1(length(Height1));            % Quota finale del fireball in km
rhof=Model_rho_atmosphere_COESA(Hf);    % Densità atmosferica alla quota finale in kg/km^3
vf=fit_velocity(length(Height1));       % Stima velocità del fireball nel punto finale in km/s (dal modello di Ceplecha o quello Lineare)
Af=final_acceleration_estimate;         % Stima accelerazione meteoroide nel punto finale del fireball in km/s^2. 
                                        % L'accelerazione deve essere minore di 0 per ottenere un modello dinamico realistico.

disp('   ')
disp(' RAW ESTIMATE STARTING PHYSICAL PARAMETERS FOR METEOROID MODEL  ')
disp('   ')

% Parametri iniziali di guess per il modello del meteoroide

% Stima di sigma, coefficiente di ablazione (s/km)^2, usando il modello di
% Ceplecha
[AminCep,II] = min(A2);                                                     % Accelerazione minima km/s^2
VminCep=fit_velocity(II);                                                   % Velocità nel punto di accelerazione minima (km/s)

if AminCep >= 0 % Nel caso l'accelerazione non stimabile con certezza
s=0.01;       % Stima iniziale fissa di Sigma, coefficiente di ablazione (s/km)^2
end

if AminCep < 0  % Nel caso l'accelerazione sia stimabile con certezza
s=(6/(VminCep*VminCep))*(1+((VminCep*VminCep)*z)/(2*Hs_best*(AminCep)));    % Coefficiente cinematico di ablazione del meteoroide (s^2/km^2)
end

Vinf=V_infty;                               % Stima della velocità all'infinito del meteoroide secondo il modello di Ceplecha o Lineare (km/s) 
Dfin=-(G*rhof*vf*vf)/(Af);                  % Rapporto M/S stimato nel punto finale del fireball in kg/km^3
D0=Dfin*exp(1)^((s/6)*(Vinf*Vinf-vf*vf));   % Stima del rapporto Massa/Diametro all'infinito (kg/km^2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENZIONE: verifica dell'esistenza di tempi di misura uguali fra due stazioni
% diverse. Se accade che dt=0 dà errore nella successiva paramfun che risolve le
% equazioni del moto con i parametri stimati del meteoroide.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_tempo=abs(diff(Tempo1)); % NOTA: diff(X) calculates differences between adjacent elements

TempoUguale = min(diff_tempo);

if TempoUguale == 0
    
    % Conteggio dei tempi consecutivi uguali e correzione aggiungendo 1/500 s al
    % secondo tempo
    count_zero=0;
    for i=1:length(diff_tempo)
        if diff_tempo(i)==0
            Tempo1(i+1)=Tempo1(i)+1/5000;
            count_zero=count_zero+1;
        end
    end
    
     disp('   ')
     disp(' WARNING FROM METEOROID DYNAMIC FUNCTION: TWO OR MORE TIMES FROM DIFFERENT STATIONS ARE IDENTICAL ')
     disp(' (But the entries in tspan must strictly increase or decrease).');
     disp('   ')
     disp([' NUMBERS OF EQUAL TIMES: ' num2str(count_zero) ]);
     disp(' CORRECTION OF THE PROBLEM BY ADDING 1/5000 s TO THE EQUAL TIMES');
     disp('   ')
   
end  

disp('   ')
disp(' COMPUTE METEOROID DYNAMIC MODEL WITH GUESS PARAMETERS  ')
disp('   ')

% Soluzione delle equazioni del moto del meteoroide con i parametri stimati
% X è una matrice composta da tre vettori colonna. Il primo riporta la velocità del fireball 
% ottenuta dal modello in funzione del tempo, la seconda la densità dell'aria e il terzo 
% la quota del fireball.

param0=[G, D0, s, Vinf, z, Hs_best, V0, rho0, H0];
X = paramfun(param0, Tempo1);
                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Ottimizzazione dei parametri stimati del fireball               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' OPTIMIZATION OF THE PHYSICAL PARAMETERS FOR METEOROID DYNAMIC MODEL  ')
disp('   ')

% Valori osservati di velocità, densità dell'aria e quote
soln(:,1)=Velocity1;                             % Velocità osservata del fireball in funzione del tempo
soln(:,2)=Model_rho_atmosphere_COESA(Height1);   % Densità dell'aria in funzione della quota
soln(:,3)=Height1;                               % Altezza osservata del fireball in funzione del tempo

% Valori inziziali dei parametri e delle condizioni di velocità, densità e
% quota.
init_cond=zeros(1,9);              % Pre-allocazione del vettore (velocizza i calcoli)
init_cond(1)=G;                    % Gamma, coefficiente di drag
init_cond(2)=D0;                   % M/S, rapporto massa - sezione meteoroide all'infinito, kg/km^2
init_cond(3)=s;                    % sigma, coefficiente di abalzione (s/km)^2
init_cond(4)=Vinf;                 % Velocità all'infinito del meteoroide km/s 
init_cond(5)=z;                    % Coseno medio della distanza zenitale della traiettoria 
init_cond(6)=Hs_best;              % Altezza di scala efficace dell'atmosfera per riprodurre il modello COESA (km)
init_cond(7)=V0;                   % Condizioni iniziali per la velocità (km/s)
init_cond(8)=rho0;                 % Condizioni iniziali per la densità dell'aria
init_cond(9)=H0;                   % Condizioni iniziali per l'altezza (km)

% Constrain sui parametri di best fit. 
% NOTA 1: G, z, Hs_best, velocità, densità e quota iniziali sono mantenuti costanti. 
% Si fa variare solo D0, s, Vinf e V0 ossia i parametri che caratterizzano il meteoroide.
%
% NOTA 2: l'ampiezza di variazione di Vinf è bene sia minore di quella di V0, in
% questo modo si è ragionevolmente sicuri che si otterrà Vinf > V0 (come è naturale attendersi dal punto di vista fisico), 
% specie se si analizza anche il tratto iniziale della traiettoria in
% modalità "delete".

% Limite inferiore dei parametri di best fit 
lb = [G, D0/5, s/5, Vinf-1.0, z, Hs_best, V0, rho0, H0];          

% Limite superiore dei parametri di best fit
ub = [G, 5*D0, 5*s, Vinf+1.0, z, Hs_best, V0, rho0, H0];      

% Numero di iterazioni (default: numero di init_cond x 100 = 900)
options = optimoptions('lsqcurvefit','MaxFunEvals',2000);                                                 

% Ottimizzazione parametri
[pbest, presnorm, presidual, exitflag, output, lambda, J] = lsqcurvefit(@paramfun, init_cond, Tempo1, soln, lb, ub, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametri di best fit del meteoroide %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=pbest(1);            % Gamma, coefficiente di drag
D0=pbest(2);           % M/S, rapporto massa - sezione meteoroide all'infinito, kg/km^2
s=pbest(3);            % sigma, coefficiente di ablazione (s/km)^2
Vinf=pbest(4);         % Velocità all'infinito del meteoroide km/s 
V_start=pbest(7);      % Velocità iniziale del meteoroide in km/s 

DI=D0*exp(1).^((s/6)*(soln(:,1).*soln(:,1)-Vinf*Vinf));   % Rapporto M/S calcolato con le velocità osservate lungo la traiettoria del fireball

if DI(1)/D0 > 1
    disp('   ')
    disp(' WARNING FROM METEOROID DYNAMIC: RELATIVE MASS/SECTION RATIO GREATHER THAN ONE!')
    disp('   ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo incertezza parametri di best fit del modello dinamico del meteoroide %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE STANDARD DEVIATION FOR BEST FIT PARAMETERS  ')
disp('   ')

DV=mean(abs(presidual(:,1)));   % Scarto media fra best fit e velocità osservate (km/s)
DH=mean(abs(presidual(:,3)));   % Scarto medio fra best fit e quote osservate (km)

% Calcolo deviazioni standard dei parametri di best fit
[sigma_G, sigma_D0, sigma_s, sigma_Vinf, sigma_V_start] = monte_carlo_best_parameters(Monte_Carlo_Number, init_cond, Tempo1, soln, lb, ub, options, DV, DH);

sigma_Vfin=mean(abs(presidual(:,1)))/sqrt(Nstaz); % Incertezza della velocità finale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stima diametro e massa meteoroide all'infinito nell'ipotesi sia una  %
% condrite                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_infinito=2000*3*D0/(4*3.5*10^12);          % Diametro del meteoroide all'infinito (m), supposto che sia una condrite con rho=3500 kg/m^3.
Massa_inf=(4*pi*3500*(D_infinito/2)^3)/3;    % Stima della massa all'infinito del meteoroide (kg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcolo quota, velocità, accelerazione, M/S nel punto finale osservato della traiettoria     %
%        usando il modello dinamico del meteoroide con i parametri ottimizzati                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE DYNAMIC METEOROID MODEL WITH BEST FIT PARAMETERS  ')
disp('   ')

% Ciclo if per evitare che l'ultimo valore di Tempo1 e Tempo siano uguali. 
% Se sono uguali allora t1=[Tempo1, Tempo(length(Tempo))] ha due elementi uguali e paramfun si blocca 
% perché il vettore temporale deve essere strettamente crescente.

if Tempo1(length(Tempo1))==Tempo(length(Tempo)) 
    
    t1=Tempo1;
    
else

t1=[Tempo1, Tempo(length(Tempo))]; % Vettore temporale usato per finalizzare il modello del meteoroide con aggiunto l'istante finale realmente osservato

end

Y0 = paramfun(pbest, t1);          % Soluzione delle equazioni del moto del meteoroide con i parametri ottimizzati

DI1=D0*exp(1).^((s/6)*(Y0(:,1).*Y0(:,1)-Vinf*Vinf));                 % Rapporto M/S calcolato con la velocità del modello dinamico lungo la traiettoria del fireball (kg/km^2)

Vfin=Y0(length(Y0(:,1)),1);                                          % Velocità nel punto finale in km/s
Hfin=Y0(length(Y0(:,3)),3);                                          % Quota nel punto finale in km
rho_fin=Model_rho_atmosphere_COESA(Hfin);                            % Densità dell'aria alla quota finale (kg/km^3)
Afin=-G*(rho_fin/D0)*Vfin*Vfin*exp(1)^((s/6)*(Vinf*Vinf-Vfin*Vfin)); % Accelerazione nel punto finale in km/s^2

M_su_S_finale=D0*exp(1)^((s/6)*(Vfin*Vfin-Vinf*Vinf));               % Rapporto Massa/sezione del meteoroide nel punto finale osservato (kg/km^2)
D_finale=2000*3*M_su_S_finale/(4*3.5*10^12);                         % Diametro del meteoroide nel punto finale (m), supposto che sia una condrite con rho=3500 kg/m^3.
Massa_finale=(4*pi*3500*(D_finale/2)^3)/3;                           % Stima della massa del meteoroide nel punto finale (kg)

% Calcolo incertezza M_su_S_finale con la propagazione degli errori
% standard

inc1=((M_su_S_finale)/D0)*sigma_D0;
inc2=((1/6)*(Vfin*Vfin-Vinf*Vinf)*(M_su_S_finale/10^6)*sigma_s);
inc3=((1/3)*s*Vfin*(M_su_S_finale/10^6)*sigma_Vfin);
inc4=((1/3)*s*Vinf*(M_su_S_finale/10^6)*sigma_Vinf);

sigma_M_su_S_finale=sqrt(inc1^2+inc2^2+inc3^2+inc4^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo del minimo dell'accelerazione dal modello dinamico del fireball %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=Model_rho_atmosphere_COESA(Y0(:,3));
Acc=-G*(rho/D0).*(Y0(:,1).*Y0(:,1)).*exp(1).^((s/6)*(Vinf*Vinf-Y0(:,1).*Y0(:,1)));

[Amin,I] = min(Acc);                                            % Accelerazione minima km/s^2
Vmin=Y0(I,1);                                                   % Velocità nel punto di accelerazione minima (km/s)
Hmin=Y0(I,3);                                                   % Quota nel punto di accelerazione minima (km)
sigma=(6/(Vmin*Vmin))*(1+((Vmin*Vmin)*z)/(2*Hs_best*(Amin)));   % Coefficiente cinematico di ablazione del meteoroide (s^2/km^2) (da confrontare con quello dinamico)
Tmin=t1(I);                                                     % Tempo in cui si ha l'accelerazione minima (s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit con il modello di Gritsevich (Solar System Research, 2007) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE PARAMETERS GRITSEVICH MODEL  ')
disp('   ')

y_obs=Height1/Hs_best;                      % Quote ridotte osservate (adimensionale)
v_obs=abs(Velocity1/Vinf);                  % Velocità ridotte osservate (adimensionale)
    
alpha0=(G*rho0_best*Hs_best)/(D0*COS_Zr);   % Stima coefficiente balistico dal modello dinamico (adimensionale)
beta0=(s*Vinf*Vinf)/2;                      % Stima parametro per la perdita di massa dal modello dinamico (adimensionale)
 
% ATTENZIONE: expint(x)=E1(x) e vale la relazione Ei(x)=-E1(-x)-i*pi
Model_Gritsevich=@(p, x)(log(p(1))+p(2)-log((-expint(-p(2))+expint(-p(2)*((x).^2)))/2));                       % Modello di Gritsevich completo

% Calcolo dei valori teorici del modello di Gritsevich
v_model=(min(v_obs):0.01:1);
y_gri_dynamics=Model_Gritsevich([alpha0 beta0], v_model);

% Massa ridotta del meteoroide in funzione della velocità ridotta (adimensionale)
% secondo il modello di Gritsevich
m_G=exp(1).^(-beta0*(1-v_model.*v_model)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ricostruzione teorica della curva di luce del bolide con magnitudine  %
% assoluta vs. tempo (supposto che sia una condrite)                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' COMPUTE THEORETICAL LIGHTCURVE  ')
disp('   ')

D_met=2000*3*DI1/(4*3.5*10^12);                                              % Diametro del meteoroide lungo la traiettoria (m), supposto che sia una condrite con rho=3500 kg/m^3.
Massa_met=(4*pi*3500*(D_met/2).^3)/3;                                        % Stima della massa del meteoroide lungo la traiettoria (kg)
DMassa_met=gradient(Massa_met, t1);                                          % Derivata temporale della massa del meteoroide (kg/s)
I_met=-0.04*(10^6)*(0.5*(Y0(:,1).^2).*DMassa_met+Massa_met.*(Y0(:,1).*Acc)); % Potenza del bolide in J/s

cost_sun=1.3608*10^3;    % Costante solare in W/m^2
mag_sun=-26.8;           % Magnitudine apparente del Sole
h_met=10^5;              % Quota di riferimento del bolide per il calcolo della magnitudine assoluta (100 km)

% Magnitudine assoluta del bolide

mag_assoluta_bolide=mag_sun-2.5*log10(I_met/(4*pi*(h_met^2)*cost_sun));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               OUTPUT DEI RISULTATI                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' OUTPUT RESULTS: SAVE IMAGES IN .BMP FORMAT AND DATA IN .TXT FILE  ')
disp('   ')

% Velocità (km/s) vs. tempo (s)                                                
string0=[fireball_name ' - Velocity vs time'];
figure

plot(Tempo, Velocity,'r.', 'MarkerSize', 16) % Valori osservati della velocità
hold on
plot(t1,Y0(:,1),'k.', 'MarkerSize', 16)    % Modello meteoroide con i parametri ottimizzati       
hold on
plot(Tempo1,X(:,1),'b--', 'LineWidth', 3)  % Modello meteoroide con i parametri stimati        
hold on
grid
legend('Observed speed','Best dynamic model', 'Guess dynamic model')
ylabel('V (km/s)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string0,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string0, 'bmp')   % Salva le figura in bmp
hold off

V_minima=min(Y0(:,1)); % Velocità minima osservata del fireball. Se V_minima > 5 km/s il meteoroide non è sopravvissuto al rientro e non c'è la fase di volo buio.

% Quota (km) vs. tempo (s)
string2=[fireball_name ' - Height vs time'];
figure
plot(Tempo, Height,'r.', 'MarkerSize', 16); % Quota osservata del meteoroide
hold on
plot(t1,Y0(:,3),'k.', 'MarkerSize', 16);   % Quota del meteoroide con i parametri ottimizzati      
hold on
plot(Tempo1,X(:,3),'b--', 'LineWidth', 3); % Quota del meteoroide con i parametri stimati        
hold on
grid
legend('Observed height','Dynamic model with best parameters', 'Dynamic model with guess parameters')
ylabel('Height (km)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string2,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string2, 'bmp')   % Salva le figura in bmp
hold off

% Velocità (km/s) vs. quota (km)
string3=[fireball_name ' - Velocity vs Height'];
figure
plot(Height, Velocity,'r.', 'MarkerSize', 16); % Curva velocità osservata
hold on
plot(Y0(:,3), Y0(:,1),'k.', 'MarkerSize', 16); % Curva velocità con i parametri ottimizzati
hold on
plot(X(:,3), X(:,1),'b--', 'LineWidth', 3); % Curva velocità con i parametri stimati
hold on
plot(Height1, fit_velocity, 'g--', 'LineWidth', 3)       % Curva delle velocità con il modello di Ceplecha (1961)
hold on
grid
legend('Observed speed','Dynamic model with best parameters', 'Dynamic model with guess parameters', [model_type ' model'])
ylabel('Velocity (km/s)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title(string3,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string3, 'bmp')   % Salva le figura in bmp
hold off

% M/S relativo vs. tempo (s)
string4=[fireball_name ' - Mass Section relative ratio vs time'];
figure
plot(t1, DI1/D0,'k.', 'MarkerSize', 16);   % Rapporto M/S teorico                     
hold on
plot(Tempo1, DI/D0,'r.', 'MarkerSize', 16);  % Rapporto M/S osservato                      
hold on
grid
legend('Best dynamic model', 'Observed M/S')
ylabel('M/S','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string4,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string4, 'bmp')   % Salva le figura in bmp
hold off

% Accelerazione vs. tempo (s)
string5=[fireball_name ' - Acceleration vs time'];
figure
plot(t1, Acc, 'k.', 'MarkerSize', 16);                        
hold on
plot(Tempo1, A2, 'r.', 'MarkerSize', 16); % Accelerazione ricavata dalle osservazioni fittando il modello di Ceplecha o lineare per la velocità
hold on
grid
legend('Best dynamic model', [model_type ' model'])
ylabel('A (km/s^2)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string5,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string5, 'bmp')   % Salva le figura in bmp
hold off

% Velocity residuals vs. tempo (s)
string6=[fireball_name ' - Velocity residuals vs time'];
figure
plot(Tempo1, presidual(:,1),'k.', 'MarkerSize', 16);                        
hold on
grid
legend('Computed from dynamical model - Observed')
ylabel('dV (km/s)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string6,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string6, 'bmp')   % Salva le figura in bmp
hold off

% Height residuals vs. tempo (s)
string7=[fireball_name ' - Height residuals vs time'];
figure
plot(Tempo1, presidual(:,3),'k.', 'MarkerSize', 16);                        
hold on
grid
legend('Computed from dynamical model - Observed')
ylabel('dH (km)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string7,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string7, 'bmp')   % Salva le figura in bmp
hold off

% Quota Gritsevich(km) vs. velocità
string8=[fireball_name ' - Gritsevich height vs velocity'];
figure
plot(v_obs, y_obs,'r.', 'MarkerSize', 16);               % Quota e velocità adimensionali osservate dal modello dinamico del meteoroide
hold on
plot(v_model, y_gri_dynamics, 'k.', 'MarkerSize', 16);   % Quota del meteoroide con il modello di Gritsevich ma alpha e beta calcolati con il modello dinamico
hold on
grid
legend('Observed values', 'Gritsevich model')
ylabel('Reduced Height (adimensional)','FontSize',20)
xlabel('Reduced Velocity (adimensional)','FontSize',20)
title(string8,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string8, 'bmp')   % Salva le figura in bmp
hold off

% Massa Gritsevich(km) vs. velocità (s)
string9=[fireball_name ' - Gritsevich mass vs velocity'];
figure
plot(v_model, m_G,'k.', 'MarkerSize', 16);   % Massa relativa del meteoroide    
hold on
grid
legend('Gritsevich model')
ylabel('Reduced Mass (adimensional)','FontSize',20)
xlabel('Reduced Velocity (adimensional)','FontSize',20)
title(string9,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string9, 'bmp')   % Salva le figura in bmp
hold off

% Magnitudine assoluta teorica del bolide vs. tempo (s)
string10=[fireball_name ' - Absolute mag vs time'];
figure
plot(t1, mag_assoluta_bolide, 'k.', 'MarkerSize', 16);                        
set(gca,'Ydir','reverse')
hold on
grid
legend('Computed from dynamical model')
ylabel('Absolute Mag (mag)','FontSize',20)
xlabel('Time (s)','FontSize',20)
title(string10,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string10, 'bmp')   % Salva le figura in bmp
hold off

% Densità aria vs. quota
string11=[fireball_name ' - Air density vs height'];
figure
plot(Height, rho_atmosphere_COESA(Height), 'r.', 'MarkerSize', 16); % Densità del modello atmosferico COESA (1976)
hold on
plot(Y0(:,3), Y0(:,2),'k.', 'MarkerSize', 16);   % Densità dell'aria dal modello dinamico con i parametri ottimizzati             
hold on
grid
legend('COESA model air density (1976)','Air density computed from dynamic model')
ylabel('Density (kg/km^3)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title(string11,'FontSize',20)
hold on
cd(working_path)
saveas(gcf, string11, 'bmp')   % Salva le figura in bmp
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Salvataggio dei risultati su file         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' SAVE THE DATA OF THE DYNAMICAL MODEL IN .TXT FILE: ' fireball_name '_Summary_Dynamics.txt'])
disp('      ');
fid = fopen([working_path '\' fireball_name '_Summary_Dynamics.txt'],'w'); % Path di salvataggio del file di dati
fprintf(fid, ['Fireball:' fireball_name] );
fprintf(fid, '    \n');
fprintf(fid, '    \n');

fprintf(fid, 'Computation mode:     %s      \n', mode);
tf = strcmp(mode,'delete');
if tf==1
fprintf(fid, 'Bin time:             0       \n');
fprintf(fid, 'Starting time (s):    %10.3f  \n', Tin);
fprintf(fid, 'Final time (s):       %10.3f  \n', Tfin);
else
fprintf(fid, 'Bin time:             %10.3f  \n', bin_time);
fprintf(fid, 'Starting time (s):    %10.3f  \n', Tin);
fprintf(fid, 'Final time (s):       %10.3f  \n', Tfin);
end

fprintf(fid, '    \n');
fprintf(fid, 'METEOROID DYNAMIC (Kalenichenko, A&A, 2006) \n');
fprintf(fid, 'Mean drag coefficient:                                                    %10.3f\n', G);
fprintf(fid, 'Mean drag coefficient uncertainty (if sigma=0, G is fixed):               %10.3f\n', sigma_G);
fprintf(fid, 'Pre-atmospheric velocity (km/s):                                          %10.3f\n', Vinf);
%fprintf(fid, 'Pre-atmospheric velocity uncertainty (km/s):                              %10.3f\n', mean(abs(presidual(:,1)))/sqrt(Nstaz));
fprintf(fid, 'Pre-atmospheric velocity uncertainty (km/s):                              %10.3f\n', sigma_Vinf); 
fprintf(fid, 'Pre-atmospheric velocity %s model (km/s):                           %10.2f\n', model_type, V_infty);
fprintf(fid, 'Pre-atmospheric velocity uncertainty %s model (km/s):               %10.2f\n', model_type, DV_infty);
fprintf(fid, 'Mean ablation factor (s^2/km^2):                                          %10.4f\n', s);
fprintf(fid, 'Mean ablation factor uncertainty (s^2/km^2):                              %10.4f\n', sigma_s);
fprintf(fid, 'Pre-atmospheric mass section ratio, M/S  (kg/m^2):                        %10.3f\n', D0/(10^6));
fprintf(fid, 'Pre-atmospheric mass section ratio uncertainty, M/S  (kg/m^2):            %10.3f\n', sigma_D0);
fprintf(fid, 'Terminal mass section ratio, M/S (kg/m^2):                                %10.3f\n', M_su_S_finale/(10^6));  
fprintf(fid, 'Terminal mass section ratio uncertainty, M/S  (kg/m^2):                   %10.3f\n', sigma_M_su_S_finale);
fprintf(fid, 'Drag x terminal mass section ratio, Gamma*(S/M) (m^2/kg):                 %10.6f\n', -(10^6)*Afin/(rho_fin*Vfin*Vfin));
fprintf(fid, '    \n');
fprintf(fid, 'START POINT - DYNAMICAL MODEL\n');
fprintf(fid, 'Height in the start point (km):                                            %10.3f\n', H0);
fprintf(fid, 'Velocity in the start point (km/s):                                        %10.3f\n', V_start);
fprintf(fid, 'Velocity in the start point uncertainty (km/s):                            %10.3f\n', sigma_V_start);
fprintf(fid, 'Air density in the start point (kg/km^3):                                  %10.3f\n', rho0);
fprintf(fid, '    \n');
fprintf(fid, 'RESIDUALS - DYNAMICAL MODEL\n');
fprintf(fid, 'Average residual of the speed fit (km/s):                                 %10.3f\n', mean(abs(presidual(:,1))));
fprintf(fid, 'Average residual of air density fit (kg/m^3):                             %10.3f\n', mean(abs(presidual(:,2)))/10^9);
fprintf(fid, 'Average residual of the heights fit (km):                                 %10.3f\n', mean(abs(presidual(:,3))));
fprintf(fid, '    \n');
fprintf(fid, 'MINIMUM ACCELERATION POINT - DYNAMICAL MODEL\n');
fprintf(fid, 'Time of minimum acceleration (s):                                         %10.2f\n', Tmin);
fprintf(fid, 'Minimum acceleration (km/s^2):                                            %10.2f\n', Amin);
fprintf(fid, 'Velocity at minimum acceleration (km/s):                                  %10.2f\n', Vmin);
fprintf(fid, 'Height at minimum acceleration (km):                                      %10.2f\n', Hmin);
fprintf(fid, 'Mean ablation factor computed with the minimum acceleration (s^2/km^2):   %10.4f\n', sigma);

if sigma < 0
fprintf(fid, 'WARNING: mean ablation factor negative! Badly velocity at minimum acceleration point. \n');
end

fprintf(fid, '    \n');
fprintf(fid, 'FINAL POINT - DYNAMICAL MODEL\n');
fprintf(fid, 'Height of the final point (km):               %10.3f\n', Hfin);
fprintf(fid, 'Uncertainty of the height (km):               %10.3f\n', mean(abs(presidual(:,3)))/sqrt(Nstaz));
fprintf(fid, 'Velocity in the final point (km/s):           %10.3f\n', Vfin);
fprintf(fid, 'Uncertainty of the velocity (km/s):           %10.3f\n', sigma_Vfin);
fprintf(fid, 'Acceleration in the final point (km/s^2):     %10.3f\n', Afin);
fprintf(fid, '    \n');
fprintf(fid, 'CRITERIA FOR NON ROTATING METEORITE-PRODUCING FIREBALL (Gritsevich, ASR, 2009)\n');
fprintf(fid, 'Log balistic coefficiente (from dynamical model), ln(alpha)                         %10.2f\n', log(alpha0));
fprintf(fid, 'Log mass loss coeffcient (from dynamical model), ln(beta)                           %10.2f\n', log(beta0));
fprintf(fid, '    \n');
fprintf(fid, 'UNDER THE HYPOTHESIS OF A CHONDRITIC METEORITE\n');
fprintf(fid, 'Starting diameter (m):                               %10.5f\n', D_infinito);
fprintf(fid, 'Starting mass (kg):                                  %10.5f\n', Massa_inf);
fprintf(fid, 'Final diameter (m):                                  %10.5f\n', D_finale);
fprintf(fid, 'Final mass (kg):                                     %10.5f\n', Massa_finale);
fprintf(fid, 'Absolute magnitude at luminosity maximum (mag):      %10.2f\n', min(mag_assoluta_bolide));
fprintf(fid, '   \n');
fprintf(fid, '   \n');       
fprintf(fid, 'BEST FIT MODEL DATA \n');
fprintf(fid, '   \n');
fprintf(fid, 'Time (s)      Height (km)        V (km/s)      Acc. (km/s^2)   M/S (kg/m^2)\n');

for j=1:length(t1)
        
fprintf(fid,'%3.3f \t\t\t %3.2f \t\t\t %3.2f \t\t\t %3.2f \t\t\t %5.1f \n', t1(j), Y0(j,3), Y0(j,1), Acc(j), DI1(j)/10^6);
        
end

% Chiusura del file di output
status = fclose(fid);

disp([' SAVE THE DATA TO COMPUTE DARK FLIGHT IN .TXT FILE: ' fireball_name '_Dynamics.txt'])
disp('      ');
fid1 = fopen([working_path '\' fireball_name '_Dynamics.txt'],'w'); % Path di salvataggio del file di dati
fprintf(fid1, '%% Parametri di input per il calcolo del dark fligt semplificato\n');
fprintf(fid1, '%10.4f %% Latitudine punto finale della traiettoria\n', Lat(length(Lat)));                                       % Latitudine punto finale della traiettoria
fprintf(fid1, '%10.4f %% Longitudine punto finale della traiettoria\n', Long(length(Long)));                                    % Longitudine punto finale della traiettoria
fprintf(fid1, '%10.4f %% Azimut (N -> E) del radiante apparente del meteoroide nell''ultimo punto osservato\n', C*Azimut_fin);  % Azimut (N -> E) del radiante apparente del meteoroide nell'ultimo punto osservato
fprintf(fid1, '%10.4f %% Distanza zenitale media della traiettoria (gradi)\n', C*acos(z));                                      % Distanza zenitale media della traiettoria (gradi)
fprintf(fid1, '%10.4f %% Velocità nel punto finale (km/s)\n', Vfin);                                                            % Velocità nel punto finale (km/s)
fprintf(fid1, '%10.4f %% Quota del punto finale (km)\n', Hfin);                                                                 % Quota del punto finale (km)
fprintf(fid1, '%10.4f %% Accelerazione del punto finale (km/s^2) \n', Afin);                                                    % Accelerazione del punto finale (km/s^2)
fprintf(fid1, '%10.4f %% Coefficiente di drag fissato (adimensionale) \n', G);                                                  % Coefficiente di drag (adimensionale)
fprintf(fid1, '%10.4f %% Rapporto finale M/S massa/sezione (kg/m^2)  \n', M_su_S_finale/(10^6));                                % Rapporto finale massa/sezione (kg/m^2) 
fprintf(fid1, '%10.6f %% Prodotto Drag*(S/M) finale (m^2/kg)   \n', -(10^6)*Afin/(rho_fin*Vfin*Vfin));                          % Prodotto Drag*(S/M) finale (m^2/kg)
status = fclose(fid1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dati di output per l'orbita eliocentrica %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigmaVinf = mean(abs(presidual(:,1)))/sqrt(Nstaz); % Incertezza media della velocità all'infinito in km/s dai residui del modello dinamico

sigmaVinf = sigma_Vinf;                              % Incertezza della velocità all'infinito in km/s dal Monte Carlo

sigmaHfin=mean(abs(presidual(:,3)))/sqrt(Nstaz);     % Incertezza media della quota finale in km

disp('   ')
disp(' END COMPUTATION DYNAMICAL MODEL  ')
disp('   ')

end

 

