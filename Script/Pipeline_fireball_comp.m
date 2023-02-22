% =========================================================================
% PIPELINE PER L'ELABORAZIONE AUTOMATICA DELLE OSSERVAZIONI DEI FIREBALL
%                MuFiS = Multipurpose Fireball Software
% =========================================================================
%
%                 ****** VERSIONE PER LA COMPILAZIONE ******
%
% Dopo avere chiesto i dati di input essenziali (fra cui il formato dei dati, PRISMA oppure standard), 
% lancia gli script di elaborazione di un bolide in sequenza, passando i dati da uno script a quello successivo: 
% traiettoria, modello dinamico, dark flight e orbita. 
%
% Ci sono due modalità di utilizzo: 
%   
% 1) BASIC: si usa per una prima analisi del bolide, quando si hanno osservazioni 
%    visuali o quando i dati sulla velocità sono rumorosi. 
%
% 2) ADVANCED: si usa quando le osservazioni sono digitali e a bassa incertezza 
%    oppure quando si vuole un modello fisico del meteoroide (massa, dimensioni, ablazione)
% 	 e la delimitazione al suolo dell'eventuale strewn field. 
%
% NOTA: per entrambe le modalità l'orbita eliocentrica viene calcolata solo se la velocità pre-atmosferica è superiore 
%       agli 11 km/s, mentre il volo buio viene stimato solo se la velocità finale è compresa fra 1 e 6 km/s. 
% 	  
% La modalità BASIC è quella di default. Si può ottenere la triangolazione completa della traiettoria 
% (quota e velocità in funzione del tempo, coordinate del radiante apparente ecc), il fit 
% della velocità vs. tempo con il modello cinematico di Ceplecha (1961),
% una prima ricostruzione dell'orbita eliocentrica con stima dell'incertezza fatte mediante 
% Monte Carlo e, se la velocità cinematica finale ottenuta con il modello di Ceplecha è al di 
% sotto dei 6 km/s (ma superiore a 1 km/s), un rozzo modello cinematico 
% del dark flight con relativa stima del punto di caduta al suolo dell'eventuale meteoroide residuo. 
% Come modello atmosferico per il dark flight in modalità basic si usa la legge isoterma delle atmosfere, 
% ossia un esponenziale con altezza di scala di 7.64 km. Con i soli dati cinematici la posiizone
% dello strewn field è molto incerta, anche di un fattore 2, quindi non è consigliato andare alla 
% ricerca di eventuali meteoriti sulla sola base della analisi cinematica. 
% L'orbita eliocentrica viene calcolata solo se la velocità all'infinito
% pre-atmosferica è superiore agli 11 km/s
% 	  
% Nella modalità ADVANCED MuFiS esegue sempre la triangolazione come nella modalità BASIC, ma a questa 
% fa seguire l'analisi dinamica del moto del bolide in atmosfera per ricavare velocità 
% all'infinito, massa, dimensioni, coefficiente di abalzione ecc. Come modello atmosferico per il 
% modello dinamico si usa il COESA (1976). Il coefficiente di drag è fissato dall'utente, 
% di default è 0.58. Se la velocità finale è al di sotto dei 6 km/s (ma superiore a 1 km/s),
% fa un'analisi di base del dark flight stimando la zona di caduta del probabile meteoroide 
% residuo. Nel modello del dark flight i dati di input vengono dal modello dinamico 
% e non più da quello cinematico del bolide, quindi lo strewn field è più accurato. Tuttavia 
% nel dark flight non si tiene conto dello stato effettivo dell'atmosfera, come densità, pressione, 
% temperatura e venti anche se si usa il modello COESA del 1976 (e non più il modello isotermico esponenziale). 
% Lo strewn field che ne esce è solo indicativo. Per calcoli più accurati PRISMA usa un altro software 
% in grado di modellizzare meglio la fase di dark flight. Viene calcolata anche l'orbita eliocentrica, ma usando 
% la velocità all'infinito del modello dinamico. L'orbita eliocentrica viene calcolata solo se la velocità all'infinito
% pre-atmosferica è superiore agli 11 km/s
%
% Se non si sceglie nessuna delle due modalità previste MuFiS termina la sua
% esecuzione.
%
% ATTENZIONE1: nella modalità BASIC le cartelle per dinamica e volo buio 
% non vengono create perché sarebbero inutili (le osservazioni visuali sono 
% troppo imprecise e la pipeline si bloccherebbe) ed è uno spreco di spazio. 
%
% ATTENZIONE2: se nella dinamica atmosferica l'accelerazione stimata nel punto finale della 
% traiettoria del fireball è nulla oppure positiva la pipeline esce dall'esecuzione e non 
% va avanti con il modello dinamico atmosferico, dark flight e orbita. Questo per evitare che la 
% non-convergenza del modello dinamico blocchi tutta la pipeline.
%
%
%                    Albino Carbognani (OAVdA)
%                    Versione 28 aprile 2020


clear all
opengl software % That will tell MATLAB to stop trying to use the graphics card, and use a software version of the OpenGL library
set(0,'DefaultFigureVisible','on')

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                                                                     %')
disp('%      MuFiS - Multipurpose Fireball Software  - PRISMA PROJECT       %')
disp('%                                                                     %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                                                                     %')
disp('%                     by Albino Carbognani (INAF-OAS)                 %')
disp('%                                Ver 3.00.00                          %')
disp('%                                                                     %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('   ')
disp('   ')

% ************************************************************************
% ************************** INPUT DEI DATI ******************************
% ************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Input del nome della cartella con la sigla identificativa del fireball %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
prompt = ' What is the name of the fireball folder (e.g. IT20180822)? ';
disp('   ')
fireball_name = input(prompt,'s');
disp('   ')
 
 % Default behavior if the user input is empty
 if isempty(fireball_name)
    disp('   ')
    disp(' Exit MuFiS: no name of the fireball folder! ')
    disp('   ')
    exit
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scelta del formato per i dati di input %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input_format = tipo di input, 'P' (PRISMA) = nuovo formato, 'S' (standard) =
% vecchio formato, valido fino a IT20190222.
 
disp('   ')
prompt = 'Input data format (PRISMA or standard, P/S)? ';
disp('   ')
input_format = input(prompt,'s');
disp('   ')

% Default behavior if the user input is empty
 if isempty(input_format)
    disp('   ')
    disp(' Exit MuFiS: no format for the fireball data! ')
    disp('   ')
    exit
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scelta della modalità di elaborazione dei dati %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp(' INFO:   ')
disp('   ')
disp(' BASIC: triangulation, orbit and dark flight (if applicable).')
disp(' Orbit and dark flight from kinematics.')
disp('   ')
disp(' ADVANCED: triangulation, dynamical model, simple dark flight (if applicable) and orbit. ')
disp(' Orbit and dark flight from dynamical model.')
disp('   ')
prompt = ' Basic computations? Y/N (default: Y) ';
disp('   ')
basic = input(prompt,'s');

% Default behavior if the user input is empty
 if isempty(basic)
    basic = 'Y';
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Scelta del data processing mode nel modello dinamico del fireball.
%
%   1-mode='delete', se si vogliono cancellare dei punti all'inizio o alla fine della 
%     traiettoria osservata.
%   2-mode='smooth', se si vuole ottenere una traiettoria media binnando i
%     dati ad intervalli di 0.1 s. 
% 
% In questa seconda modalità il software calcola anche l'accelerazione osservata in funzione 
% del tempo e mostra/salva le relative figure.
%                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTA: questi dati servono solo se non si tratta di osservazioni visuali! 

if (basic == 'N') || (basic == 'n')

disp('   ')
prompt = ' Data processing mode in the meteoroid_dynamic function: delete or smooth (default delete)? ';
disp('   ')
mode = input(prompt, 's');
tf = strcmp(mode, 'smooth'); % Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.

    if tf==1
        
       mode='smooth'; start_delay = 0;  end_delay = 0;
       disp('   ')
       disp( ' Data processing in smooth mode. Starting/Ending delay: 0 ');
       disp('   ')
       prompt = ' Bin time to avoid noisy data (default: 0.075 s)? ';
       disp('   ')
       bin_time = input(prompt);
           
       % Default behavior if the user input is empty
       if isempty(bin_time)
          bin_time=0.075;
       end
       
    else
           mode='delete';
           disp('   ')
           prompt = ' Starting delay to avoid noisy data (fraction of the total flight time, default: 0.05)? ';
           disp('   ')
           start_delay = input(prompt);

           % Default behavior if the user input is empty
           if isempty(start_delay)
               start_delay=0.05;
           end

           disp('   ')
           prompt = ' Ending delay to avoid noisy data (fraction of the total flight time, default: 0.05)? ';
           disp('   ')
           end_delay = input(prompt);

           % Default behavior if the user input is empty
           if isempty(end_delay)
              end_delay=0.05;
           end
       
           bin_time=0;
    end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scelta del coefficiente di drag %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('   ')
    prompt = ' Fixed drag coefficient for the fireball phase (default: 0.58)? ';
    disp('   ')
    G = input(prompt);   
     
    % Default behavior if the user input is empty
       if isempty(G)
          G=0.58;
       end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting dei paths di lavoro %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
prompt = ' Full path of the fireball folder (e.g. C:/Users/Name/Documents/MATLAB/Fireball_data/)? ';
disp('   ')
data_path = input(prompt,'s');

% Default behavior if the user input is empty
 if isempty(data_path)
    disp('   ')
    disp(' Exit MuFiS: no data path! ')
    disp('   ')
    exit 
 end

%%% data_path='C:\Users\demo\Documents\MATLAB\Fireball\Dati_bolidi'; % Path per la cartella con i dati di tutti i bolidi
path1=[data_path '\' fireball_name '\Data'];                         % Cartella con i dati osservati del fireball

%%% script_path='C:\Users\demo\Documents\MATLAB\Fireball';           % Path per la cartella contenente gli script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creazione delle cartelle per ospitare i risultati sulla dinamica del meteoroide, dark flight e orbita 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' INFO: Creation of the folders to host results about trajectory, dynamics, dark flight and orbit. ')
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
%%% addpath([script_path '\Multiple_Trajectories'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
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
    % dell'orbita eliocentrica viene fatto solo se la velocità all'infinito del meteoroide è
    % superiore a 11 km/s.
       
    if (V_infty > 11) % Velocità all'infinito del fireball ottenuta dal modello cinematico di Ceplecha in km/s
    
    %addpath([script_path '\Orbit'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
    Orbit_pipeline(path5, fireball_name, jd0, UT, AR_Rad, Dec_Rad, inc_rad, Lat_fin, Long_fin, V_infty, DV_infty, quota_finale, d_quota_finale, 100); 
    
    disp('   ')
    disp(' BASIC COMPUTATION: EXIT MuFiS AFTER TRIANGULATION AND ORBIT, NO DARK-FLIGHT ')
    disp('   ')
    prompt = ' Press Enter key to exit '; % In questo modo le figure restano aperte alla fine dei calcoli
    input(prompt);
    disp('   ')
    exit
    
    else
        
    disp('   ')
    disp(' WARNING BASIC COMPUTATION: Meteoroid to slow, i.e. V_infty < 11 km/s, to compute heliocentric orbit!')
    disp('   ')    
    prompt = ' Press Enter key to exit ';
    input(prompt);
    disp('   ')
    exit  

    end
    
end

if ((basic == 'Y') || (basic == 'y')) && (V_finale <= 6.0) && (V_finale > 1.0)
    
    % Calcolo dell'orbita e del dark flight con i soli dati della triangolazione. 
    % In questo modo si ha un'idea dell'orbita e di dove sia finito il meteoroide residuo senza dover passare 
    % attraverso il calcolo del modello dinamico del bolide. Le incertezze orbitali sono date direttamente dalla triangolazione invece che dal modello
    % dinamico. Vengono calcolati 100 cloni Monte Carlo. Il calcolo dell'orbita eliocentrica viene fatto solo se la velocità all'infinito del meteoroide è
    % superiore a 11 km/s.
    
    if (V_infty > 11) % Velocità all'infinito del fireball ottenuta dal modello cinematico di Ceplecha in km/s
    
      %addpath([script_path '\Orbit'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
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
        
    %addpath([script_path '\Dark_flight'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
    Simple_Dark_flight_meteoroid_basic(path4, fireball_name, Lat_fin, Long_fin, Azimut_finale, 90-inclinazione_finale, V_finale, quota_finale, A_finale);
    
    disp('   ')
    disp(' BASIC COMPUTATION: EXIT MuFiS AFTER TRIANGULATION, ORBIT AND DARK-FLIGHT PATH ESTIMATE ')
    disp('   ')
    prompt = ' Press Enter key to exit '; % In questo modo le figure restano aperte alla fine dei calcoli
    input(prompt);
    disp('   ')
    exit
    
end

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

%%% addpath([script_path '\Meteoroid_dynamic'])                                         % Adds the specified folders to the top of the search path for the current MATLAB® session.
Tin=start_delay*fireball_total_duration; Tfin=(1-end_delay)*fireball_total_duration;    % Starting/ending time in s for fireball dynamic analysis
                                                                                        % Nstaz=stations number
                                                                        
disp('   ')
prompt = 'Number of trajectory scenarios to compute with Monte Carlo method (default = 100)? ';
disp('   ')
Nmc = input(prompt);
 
% Default parameter if the user input is empty
if isempty(Nmc)
    Nmc=100; % Numero di scenari da calcolare con il metodo Monte Carlo
end                                                                                        
                                                                                        
[V_inf, SV_inf, quota_fin, Squota_fin, V_minima, acc_fin_est]=Meteoroid_dynamic(path2, path3, fireball_name, Nstaz, Tin, Tfin, Azimut, mode, G, bin_time, Nmc); 

% La pipeline termina se il valore stimato della accelerazione del bolide
% nel punto finale è positiva. Questo per evitare la non convergenza del
% modello dinamico con conseguente errore.

if acc_fin_est >= 0
    disp('   ')
    disp(' Exit MuFiS. Possible reasons: ')
    disp('   ')
    disp(' 1-Raw estimate of the final acceleration in "Meteoroid Dynamic" is zero or positive! ')
    disp('   ')
    disp(' 2-Two observed time from different stations are identical but vector time in "Meteoroid Dynamic" must strictly increase! ')
    disp('   ')
    prompt = ' Press Enter key to exit ';
    input(prompt);
    disp('   ')
    exit 
end

% ATTENZIONE: la fase di volo buio e il punto di impatto al suolo vengono calcolati solo se la velocità minima in atmosfera del fireball 
% (fornita dal modello dinamico) è compresa fra 2 e 6 km/s, ossia solo se è abbastanza probabile (tenuto conto che la velocità ha un'incertezza di circa 1-2 km/s)
% che esista un meteoroide residuo, alias meteorite al suolo.

if (V_minima <= 6) && (V_minima > 1.0) % Velocità minima del fireball ottenuta dal modello dinamico in km/s

disp('   ')
disp(' WARNING ADVANCED COMPUTATION: Meteoroid survived fireball phase (1 km/s < Vmin < 6 km/s)!')
disp('   ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dark flight e punto di impatto al suolo %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% addpath([script_path '\Dark_flight'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
Simple_Dark_flight_meteoroid(path3, path4, fireball_name); 

else
    disp('   ')
    disp(' WARNING ADVANCED COMPUTATION: Meteoroid too fast (Vmin > 6 km/s): no dark flight and strewn field computed')
    disp('   ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dell'orbita eliocentrica (avviene solo se la velocità all'infinito è pari o superiore a 11 km/s)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (V_inf >= 11) % Velocità all'infinito del fireball ottenuta dal modello dinamico in km/s

 %%% addpath([script_path '\Orbit'])  % Adds the specified folders to the top of the search path for the current MATLAB® session.
  
 disp('   ')
 prompt = ' Number of heliocentric orbits to compute with Monte Carlo method (default = 100)? ';
 disp('   ')
 N = input(prompt);
 
 % Default parameter if the user input is empty
 if isempty(N)
     N=100; % Numero di orbite di default da calcolare con il metodo Monte Carlo
 end

 Orbit_pipeline(path5, fireball_name, jd0, UT, AR_Rad, Dec_Rad, inc_rad, Lat_fin, Long_fin, V_inf, SV_inf, quota_fin, Squota_fin, N); 

 disp('   ')
 disp( ' End advanced fireball computations. Congratulations! ');
 disp('   ')
 
 % result = input('Press Enter key to exit'); % In questo modo la finestra DOS resta aperta alla fine dei calcoli
 
 disp('   ') 
 prompt = ' Press Enter key to exit ';
 input(prompt);
 disp('   ')
 exit
 
else
   
disp('   ')
disp(' WARNING ADVANCED COMPUTATION: Meteoroid too slow, i.e. V_inf < 11 km/s, to compute heliocentric orbit!')
disp('   ')    
prompt = ' Press Enter key to exit ';
input(prompt);
disp('   ')
exit  
    
end

% ************************************** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         END ADVANCED COMPUTATION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************************** %

else
    
   disp('   ')
   disp(' No valid computation mode (basic or advanced) selected  ') 
   prompt = ' Press Enter key to exit ';
   input(prompt);
   disp('   ')
   exit
 
end
