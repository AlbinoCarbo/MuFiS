% =============================================================================
% SUPER PIPELINE PER L'ELABORAZIONE AUTOMATICA DELLE OSSERVAZIONI DEI FIREBALL
%                   MuFiS = Multipurpose Fireball Software
% =============================================================================
%
%                 ****** VERSIONE PER LA COMPILAZIONE ******
%
% Dopo avere letto i dati di input essenziali dal file "Fireballs_list.txt", lancia gli script 
% di elaborazione per N bolidi in sequenza, passando i dati da uno script a quello successivo: 
% traiettoria, modello dinamico, dark flight e orbita. In questo modo si possono elaborare grandi quantit� di bolidi.
% L'output delle figure � soppresso.
%
% Ci sono due modalit� di utilizzo: 
%   
% 1) BASIC: si usa per una prima analisi del bolide, quando si hanno osservazioni 
%    visuali o quando i dati sulla velocit� sono rumorosi. 
%
% 2) ADVANCED: si usa quando le osservazioni sono digitali e a bassa incertezza 
%    oppure quando si vuole un modello fisico del meteoroide (massa, dimensioni, ablazione)
% 	 e la delimitazione al suolo dell'eventuale strewn field. 
%
% NOTA: per entrambe le modalit� l'orbita eliocentrica viene calcolata solo se la velocit� pre-atmosferica � superiore 
%       agli 11 km/s, mentre il volo buio viene stimato solo se la velocit� finale � compresa fra 1 e 6 km/s. 
% 	  
% La modalit� BASIC � quella di default. Si pu� ottenere la triangolazione completa della traiettoria 
% (quota e velocit� in funzione del tempo, coordinate del radiante apparente ecc), il fit 
% della velocit� vs. tempo con il modello cinematico di Ceplecha (1961),
% una prima ricostruzione dell'orbita eliocentrica con stima dell'incertezza fatte mediante 
% Monte Carlo e, se la velocit� cinematica finale ottenuta con il modello di Ceplecha � al di 
% sotto dei 6 km/s (ma superiore a 1 km/s), un rozzo modello cinematico 
% del dark flight con relativa stima del punto di caduta al suolo dell'eventuale meteoroide residuo. 
% Come modello atmosferico per il dark flight in modalit� basic si usa la legge isoterma delle atmosfere, 
% ossia un esponenziale con altezza di scala di 7.64 km. Con i soli dati cinematici la posiizone
% dello strewn field � molto incerta, anche di un fattore 2, quindi non � consigliato andare alla 
% ricerca di eventuali meteoriti sulla sola base della analisi cinematica. 
% L'orbita eliocentrica viene calcolata solo se la velocit� all'infinito
% pre-atmosferica � superiore agli 11 km/s
% 	  
% Nella modalit� ADVANCED MuFiS esegue sempre la triangolazione come nella modalit� BASIC, ma a questa 
% fa seguire l'analisi dinamica del moto del bolide in atmosfera per ricavare velocit� 
% all'infinito, massa, dimensioni, coefficiente di abalzione ecc. Come modello atmosferico per il 
% modello dinamico si usa il COESA (1976). Il coefficiente di drag � fissato dall'utente, 
% di default � 0.58. Se la velocit� finale � al di sotto dei 6 km/s (ma superiore a 1 km/s),
% fa un'analisi di base del dark flight stimando la zona di caduta del probabile meteoroide 
% residuo. Nel modello del dark flight i dati di input vengono dal modello dinamico 
% e non pi� da quello cinematico del bolide, quindi lo strewn field � pi� accurato. Tuttavia 
% nel dark flight non si tiene conto dello stato effettivo dell'atmosfera, come densit�, pressione, 
% temperatura e venti anche se si usa il modello COESA del 1976 (e non pi� il modello isotermico esponenziale). 
% Lo strewn field che ne esce � solo indicativo. Per calcoli pi� accurati PRISMA usa un altro software 
% in grado di modellizzare meglio la fase di dark flight. Viene calcolata anche l'orbita eliocentrica, ma usando 
% la velocit� all'infinito del modello dinamico. L'orbita eliocentrica viene calcolata solo se la velocit� all'infinito
% pre-atmosferica � superiore agli 11 km/s
%
% ATTENZIONE1: nella modalit� BASIC le cartelle per dinamica e volo buio 
% non vengono create perch� sarebbero inutili (le osservazioni visuali sono 
% troppo imprecise e la pipeline si bloccherebbe) ed � uno spreco di spazio. 
%
% ATTENZIONE2: se nella dinamica atmosferica l'accelerazione stimata nel punto finale della 
% traiettoria del fireball � nulla oppure positiva la pipeline esce dall'esecuzione e non 
% va avanti con il modello dinamico atmosferico, dark flight e orbita. Questo per evitare che la 
% non-convergenza del modello dinamico blocchi tutta la pipeline.
%
%
%                    Albino Carbognani (INAF-OAS)
%                    Versione 28 aprile 2020


clear all
opengl software % That will tell MATLAB to stop trying to use the graphics card, and use a software version of the OpenGL library

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                                                                     %')
disp('%      MuFiS - Multipurpose Fireball Software  - PRISMA PROJECT       %')
disp('%                                                                     %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                       PIPELINE VERSION                              %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                                                                     %')
disp('%               by Albino Carbognani (INAF-OATo/OAVdA)                %')
disp('%                             Ver 1.10.00                             %')
disp('%                                                                     %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('   ')
disp('   ')

% ************************************************************************
% ************************** INPUT DEI DATI ******************************
% ************************************************************************

% Soppressione della visualizzazione delle figure
set(0,'DefaultFigureVisible','off')

% Load tabella con l'elenco dei fireball da elaborare
s0='Fireballs_list.txt';
T = readtable(s0);

% Inizializzazione array di stringhe o numeri
FireballName=table2array(T(:,1));       % Array di stringhe
Format=table2array(T(:,2));             % Array di stringhe
Basic=table2array(T(:,3));              % Array di stringhe
Mode=table2array(T(:,4));               % Array di stringhe
BinTime=table2array(T(:,5));            % Array numerico
Drag=table2array(T(:,6));               % Array numerico
Tinizio=table2array(T(:,7));            % Array numerico
Tfine=table2array(T(:,8));              % Array numerico
MonteCarloOrbit=table2array(T(:,9));    % Array numerico
Nmc=table2array(T(:,10));               % Array numerico

NumeroFireballs=length(FireballName);

% Setting del path dove risiedono le cartelle con i dati dei bolidi

disp('   ')
prompt = ' Full path of the fireballs folders (e.g. C:/Users/Name/Documents/MATLAB/Fireball_data/)? ';
disp('   ')
data_path = input(prompt,'s');

% Default behavior if the user input is empty
 if isempty(data_path)
    disp('   ')
    disp(' Exit MuFiS: no data path! ')
    disp('   ')
    exit 
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inizio del ciclo di elaborazione per tutti i fireball osservati %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:NumeroFireballs
    
% Inizializzazione variabili per il singolo bolide    
    
fireball_name = char(FireballName(ii));   % Nome dell'ii-esimo fireball

input_format = char(Format(ii));          % Formato dei dati per l'ii-esimo fireball

basic = char(Basic(ii));                  % Modalit� di elaborazione dell'ii-esimo fireball

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Scelta del data processing mode nel modello dinamico del fireball.
%
%   1-mode='delete', se si vogliono cancellare dei punti all'inizio o alla fine della 
%     traiettoria osservata.
%   2-mode='smooth', se si vuole ottenere una traiettoria media binnando i
%     dati ad intervalli definiti. 
% 
% In questa seconda modalit� il software calcola anche l'accelerazione osservata in funzione 
% del tempo e mostra/salva le relative figure.
%                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = char(Mode(ii));                     % Setting dinamica atmosferica dell'ii-esimo fireball

bin_time=BinTime(ii);                      % Ampiezza del bin temporale dell'ii-esimo fireball
  
G=Drag(ii);                                % Coefficiente di drag dell'ii-esimo fireball

start_delay=Tinizio(ii);                   % Frazione del tempo da scartare dopo l'inizio dell'ii-esimo fireball

end_delay=Tfine(ii);                       % Frazione del tempo da scartare prima della fine dell'ii-esimo fireball

NumberMonteCarloOrbit=MonteCarloOrbit(ii); % Numero di orbite dei cloni Monte Carlo

disp('   ')
disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('   ')
disp([' STARTING COMPUTATION FOR FIREBALL: ' fireball_name ])
disp('   ') 
disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('   ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting del path di lavoro  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path1=[data_path '\' fireball_name '\Data'];                         % Cartella con i dati osservati del fireball

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creazione delle cartelle per ospitare i risultati sulla dinamica del meteoroide, dark flight e orbita 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp([' INFO ' fireball_name ': Creation of the folders to host results about trajectory, dynamics, dark flight and orbit. '])
disp(' They are created inside fireball folder and if they don''t exist only. ')
disp('   ')

cd([data_path '\' fireball_name])
  
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir 'Trajectory'
mkdir 'Orbit';

if (basic == 'N') || (basic == 'n')
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  mkdir 'Dynamic'
  mkdir 'Dark_flight';
  mkdir 'Orbit';
end
  
% Definizione dei path di salvataggio risultati per traiettoria, dinamica, dark flight e orbita
path2=[data_path '\' fireball_name '\Trajectory'];
path3=[data_path '\' fireball_name '\Dynamic'];
path4=[data_path '\' fireball_name '\Dark_flight'];
path5=[data_path '\' fireball_name '\Orbit'];
  
% ************************************************************************
% ************************* FINE INPUT DEI DATI **************************
% ************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Triangolazione da N stazioni %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% Funzione accessoria: legge i dati temporali del fireball, ossia i JD, dalla stazione n.1 e calcola l'UT medio e il
% giorno giuliano per le 0 UT
[jd0, UT] = JD0_UT_function(path1, fireball_name, input_format);

% Funzione per la triangolazione da N stazioni contemporaneamente.
%%% addpath([script_path '\Multiple_Trajectories'])  % Adds the specified folders to the top of the search path for the current MATLAB� session.
[Azimut, fireball_total_duration, Nstaz, AR_Rad, Dec_Rad, Lat_fin, Long_fin, inc_rad, V_infty, DV_infty, quota_finale, d_quota_finale, V_finale, Azimut_finale, inclinazione_finale, A_finale]=Trajectory_multi_stations(path1, path2, fireball_name, input_format, jd0, UT);
    
% ************************************** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            BASIC COMPUTATION           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************************** %

if ((basic == 'Y') || (basic == 'y')) && (V_finale > 6.0)
    
    % Calcolo dell'orbita con i soli dati della triangolazione. 
    % In questo modo si ha un'idea dell'orbita eliocentrica del meteoroide progenitore senza dover passare 
    % attraverso il calcolo del modello dinamico del bolide.
    % Le incertezze sono date dalla triangolazione invece che dal modello
    % dinamico. Vengono calcolati 100 cloni Monte Carlo. Il calcolo
    % dell'orbita eliocentrica viene fatto solo se la velocit� all'infinito del meteoroide �
    % superiore a 11 km/s.
       
    if (V_infty > 11) % Velocit� all'infinito del fireball ottenuta dal modello cinematico di Ceplecha in km/s
    
    %addpath([script_path '\Orbit'])  % Adds the specified folders to the top of the search path for the current MATLAB� session.
    Orbit_pipeline(path5, fireball_name, jd0, UT, AR_Rad, Dec_Rad, inc_rad, Lat_fin, Long_fin, V_infty, DV_infty, quota_finale, d_quota_finale, 100); 
    
    disp('   ')
    disp(' BASIC COMPUTATION: TRIANGULATION AND ORBIT, NO DARK-FLIGHT ')
    disp('   ')
        
    else
        
    disp('   ')
    disp(' WARNING BASIC COMPUTATION: Meteoroid to slow, i.e. V_infty < 11 km/s, to compute heliocentric orbit!')
    disp('   ')    
    disp(' BASIC COMPUTATION: TRIANGULATION ONLY, NO DARK-FLIGHT ') 
    disp('   ')
    
    end
    
end % Fine ciclo if basic computation senza dark-flight

if ((basic == 'Y') || (basic == 'y')) && (V_finale <= 6.0) && (V_finale > 1.0)
    
    % Calcolo dell'orbita e del dark flight con i soli dati della triangolazione. 
    % In questo modo si ha un'idea dell'orbita e di dove sia finito il meteoroide residuo senza dover passare 
    % attraverso il calcolo del modello dinamico del bolide. Le incertezze orbitali sono date direttamente dalla triangolazione invece che dal modello
    % dinamico. Vengono calcolati 100 cloni Monte Carlo. Il calcolo dell'orbita eliocentrica viene fatto solo se la velocit� all'infinito del meteoroide �
    % superiore a 11 km/s.
    
    if (V_infty > 11) % V_inf=velocit� all'infinito del fireball ottenuta dal modello cinematico di Ceplecha in km/s
    
      %addpath([script_path '\Orbit'])  % Adds the specified folders to the top of the search path for the current MATLAB� session.
      Orbit_pipeline(path5, fireball_name, jd0, UT, AR_Rad, Dec_Rad, inc_rad, Lat_fin, Long_fin, V_infty, DV_infty, quota_finale, d_quota_finale, 100); 
    
    else
        
      disp('   ')
      disp(' WARNING BASIC COMPUTATION: Meteoroid too slow, i.e. V_infty < 11 km/s, to compute heliocentric orbit!')
      disp('   ') 
    
    end
    
    disp('   ')
    disp('WARNING BASIC COMPUTATION: final velocity under 6 km/s, probably meteoroid survived fireball phase!')
    disp('Computation with the atmospheric dynamic model is necessary.')
    disp('   ')
    
    % Stima del dark flight
     
    cd([data_path '\' fireball_name])
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir 'Dark_flight'; % Creazione della cartella per il dark flight
        
    %addpath([script_path '\Dark_flight'])  % Adds the specified folders to the top of the search path for the current MATLAB� session.
    Simple_Dark_flight_meteoroid_basic(path4, fireball_name, Lat_fin, Long_fin, Azimut_finale, 90-inclinazione_finale, V_finale, quota_finale, A_finale);
    
    disp('   ')
    disp(' BASIC COMPUTATION: TRIANGULATION, ORBIT AND DARK-FLIGHT PATH ESTIMATED ')
    disp('   ')
        
end % Fine ciclo if basic computation con dark-flight

% ************************************** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          END BASIC COMPUTATION         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************************** %

% ************************************** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        START ADVANCED COMPUTATION      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************************** %

if ((basic == 'N') || (basic == 'n'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analisi della dinamica atmosferica del meteoroide %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% addpath([script_path '\Meteoroid_dynamic'])                                         % Adds the specified folders to the top of the search path for the current MATLAB� session.
Tin=start_delay*fireball_total_duration; Tfin=(1-end_delay)*fireball_total_duration;    % Starting/ending time in s for fireball dynamic analysis
                                                                                        % Nstaz=stations number
                                                                        
[V_inf, SV_inf, quota_fin, Squota_fin, V_minima, acc_fin_est]=Meteoroid_dynamic(path2, path3, fireball_name, Nstaz, Tin, Tfin, Azimut, mode, G, bin_time, Nmc); 

% La pipeline termina se il valore stimato della accelerazione del bolide
% nel punto finale � positiva. Questo per evitare la non convergenza del
% modello dinamico con conseguente errore.

if acc_fin_est >= 0
    disp('   ')
    disp(' WARNING FROM ADVANCED COMPUTATION: EXIT FIREBALL PIPELINE. ')
    disp('   ')
    disp(' RAW ESTIMATE OF THE FINAL ACCELERATION IN "METEOROID DYNAMIC" IS ZERO OR POSITIVE! ')
    disp(' NO INTEGRATION OF THE DIFFERENTIAL MOTION EQUATIONS POSSIBLE. ')
    disp(' SUGGESTION: USE BASIC COMPUTATION. ')
    disp('   ')
    disp(' Final acceleration estimate (km/s^2) ')
    acc_fin_est
    return
end

% ATTENZIONE: la fase di volo buio e il punto di impatto al suolo vengono calcolati solo se la velocit� minima in atmosfera del fireball 
% (fornita dal modello dinamico) � compresa fra 1 e 6 km/s, ossia solo se � abbastanza probabile (tenuto conto che la velocit� ha un'incertezza di circa 1-2 km/s)
% che esista un meteoroide residuo, alias meteorite al suolo.

 if (V_minima <= 6) && (V_minima > 1.0) % V_minima=velocit� minima del fireball ottenuta dal modello dinamico in km/s

  disp('   ')
  disp(' WARNING FROM ADVANCED COMPUTATION: Meteoroid survived fireball phase (1 km/s < Vmin < 6 km/s)!')
  disp('   ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dark flight e punto di impatto al suolo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Simple_Dark_flight_meteoroid(path3, path4, fireball_name); 

  else
    disp('   ')
    disp(' WARNING FROM ADVANCED COMPUTATION: Meteoroid too fast (Vmin > 6 km/s): no dark flight and strewn field computed')
    disp('   ')
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dell'orbita eliocentrica (avviene solo se la velocit� all'infinito � pari o superiore a 11 km/s)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (V_inf >= 11) % V_inf=velocit� all'infinito del fireball ottenuta dal modello dinamico in km/s

  N = NumberMonteCarloOrbit;
 
  Orbit_pipeline(path5, fireball_name, jd0, UT, AR_Rad, Dec_Rad, inc_rad, Lat_fin, Long_fin, V_inf, SV_inf, quota_fin, Squota_fin, N); 
 
 else
   
  disp('   ')
  disp(' WARNING FROM ADVANCED COMPUTATION: Meteoroid too slow, i.e. V_inf < 11 km/s, to compute heliocentric orbit!')
  disp('   ')    
      
 end
 
end % Fine istruzione if per advanced computation

% ************************************** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         END ADVANCED COMPUTATION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************************** %

disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('                                                        ')
disp(['    END COMPUTATION FOR FIREBALL: ' fireball_name      ])
disp('                                                        ') 
disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('   ')

end % Fine del ciclo for dei fireballs

% Uscita dallo script

disp('   ')
disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('                                                        ')
disp('              END COMPUTATION PIPELINE                  ')
disp('                                                        ') 
disp(' *****************************************************  ')
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ')
disp(' *****************************************************  ')
disp('   ')

exit
