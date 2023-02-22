% Funzione per la stima preliminare della velocità all'infinito (con incertezza) 
% del meteoroide con il modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h).  
% Utile per la stima preliminare dell'orbita con i soli dati di
% triangolazione, senza passare dal modello dinamico che, se le
% osservazioni non sono accurate, potrebbe non convergere.
%
% Questa funzione fa anche una stima della accelerazione nel punto finale osservato della
% traiettoria, in modo tale che si possa fare una stima della lunghezza del
% dark flight senza passare per il modello con la dinamica atmosferica.
%
% ATTENZIONE: se il modello di Ceplecha non converge vuol dire che il
% fit non è buono, probabilmente perché il bolide è molto veloce e la velocità in 
% funzione della quota tende ad essere lineare invece che esponenziale. In
% questo caso si passa ad un fit lineare con la quota invece di uno
% esponenziale come quello di Ceplecha. La velocità all'infinito è quella
% che si ottiene estrapolando i valori osservati a 120 km di quota.
%
% INPUT:
% Height1 = quote osservate del bolide (km)
% Velocity1 = velocità osservate del bolide (km/s)
% inclinazione = inclinazione della traiettoria (radianti)
%
% OUTPUT:
% V_infty = velocità del bolide all'infinito (km/s)
% DV_infty = deviazione standard della velocità all'infinito (km/s)
% fit_velocity = vettore del fit delle velocità con il modello di Ceplecha
% Afin = stima accelerazione nel punto finale della traiettoria osservata
% (km/s^2)
% model_type = tipo di modello cinematico utilizzato, 'Ceplecha' oppure '  Linear'.
% In questo modo nel file di output della triangolazione viene citato il
% modello esatto da cui sono stati ottenuti i parametri cinematici.
%
% Albino Carbognani, versione del 30 aprile 2019


function [V_infty, DV_infty, fit_velocity, Afin, model_type]=Ceplecha_model_velocity(Height1, Velocity1, inclinazione)

% Definizione del modello di Ceplecha per la velocità in funzione della
% quota (Ceplecha, 1961).

Model_Ceplecha_V=@(p, x)(p(1)+p(2)*exp(1).^(-p(3)*x)); 

V2=median(Velocity1); % Valore guess della velocità all'infinito. 

% NOTA: La funzione mediana è poco sensibile ai punti più scatterati della velocità 
% e a quelli che si trovano verso la fine della traiettoria del bolide
% quando la decelerazione è elevata.

warning('off', 'MATLAB:rankDeficientMatrix')
warning('off', 'stats:nlinfit:ModelConstantWRTParam') % Sopprime i warning nel caso di fit scadente (tanto si vede dai plot del main)

disp('   ')
    disp(' From Ceplecha model velocity function: COMPUTE CEPLECHA VELOCITY MODEL  ')
    disp('   ')

startingVals = [V2 mean(Velocity1) 0.05];         % Valori iniziali stimati dei parametri p(1), p(2) e p(3)
[coefEsts, R1, J1, CovB1, MSE1] = nlinfit(Height1, Velocity1, Model_Ceplecha_V, startingVals);

fit_velocity=Model_Ceplecha_V(coefEsts, Height1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stima velocità all'infinito %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_infty=coefEsts(1); 

DV_infty=sqrt(CovB1(1,1)); % Incertezza velocità all'infinito

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stima dell'accelerazione nel punto finale della traiettoria osservata con il modello di Ceplecha 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Il valore dell'accelerazione nel punto finale osservato della traiettoria è dato da:
% a(t)=dV/dt=(dV/dh)*(dh/dt)=(-kCv*exp(-kh))*(Vh), dove Vh < 0 è la velocità di caduta lungo la verticale 

Ccep=coefEsts(2); % Secondo coefficiente modello di Ceplecha
Kcep=coefEsts(3); % Terzo coefficiente modello di Ceplecha

Vh=-(fit_velocity(length(Height1)))*sin(inclinazione);               % Stima della velocità di caduta verticale nel punto finale osservato (km/s)
Afin=-(Ccep*Kcep)*(exp(1).^-(Kcep*Height1(length(Height1)))).*Vh;    % Accelerazione finale nel modello di Ceplecha (km/s^2)

model_type='Ceplecha';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modello lineare fra velocità e quota, adatto per bolidi molto veloci   %
% che si disintegrano subito.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTA: il modello lineare entra in gioco se il modello della velocità di Ceplecha non converge entro le 100 iterazioni 

% Legge l'eventuale warning che arriva dal bad fit del modello di Ceplecha
[warnmsg, msgid] = lastwarn;                                

% Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
tf1=strcmp(msgid,'MATLAB:rankDeficientMatrix');
tf2=strcmp(msgid,'stats:nlinfit:ModelConstantWRTParam');
tf3=strcmp(msgid,'stats:nlinfit:IterationLimitExceeded');  
tf4=strcmp(msgid,'stats:nlinfit:IllConditionedJacobian');

if (tf1 == 1) || (tf2 == 1) || (tf3 == 1) || (tf4 == 1) 
    
    disp('   ')
    disp(' From Ceplecha model velocity function: BAD CEPLECHA MODEL, SKYP TO LINEAR VELOCITY MODEL  ')
    disp('   ')
    
   Model_Ceplecha_V2=@(p, x)(p(1)+p(2)*x);               % Modello lineare fra velocità e quota x 
   startingVals2 = [V2 0];                               % Valori iniziali stimati dei parametri p(1), p(2)
   [coefEsts2, R1, J1, CovB2, MSE1] = nlinfit(Height1, Velocity1, Model_Ceplecha_V2, startingVals2);

   fit_velocity=Model_Ceplecha_V2(coefEsts2, Height1);
   
   V_infty=Model_Ceplecha_V2(coefEsts2, 120);              % La velocità all'infinito del meteoroide è quella a 120 km di quota
   DV_infty=sqrt(CovB2(1,1));                              % Incertezza velocità all'infinito
   Vh=-(fit_velocity(length(Height1)))*sin(inclinazione);  % Stima della velocità di caduta verticale nel punto finale osservato (km/s)
   Afin=coefEsts2(2)*Vh;                                   % Accelerazione finale nel modello lineare (km/s^2)
   
   model_type='  Linear';
   
end

end
