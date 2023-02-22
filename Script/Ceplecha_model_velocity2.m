% Funzione per la stima preliminare della velocità all'infinito (con incertezza) 
% del meteoroide con il modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h).  
% Se il modello di Ceplecha non converge si passa ad un modello lineare.
% Questa funzione fa anche una stima dell'accelerazione in funzione del tempo.
%
% ATTENZIONE: se il modello di Ceplecha non converge vuol dire che il
% fit non è buono, probabilmente perché il bolide è molto veloce e la velocità in 
% funzione della quota tende ad essere lineare invece che esponenziale. In
% questo caso si passa ad un fit lineare con la quota invece di uno
% esponenziale come quello di Ceplecha. La velocità all'infinito è quella
% che si ottiene estrapolando i valori osservati a 120 km di quota.
%
% INPUT:
% Height = vettore di tutte le quote osservate del bolide (km)
% Height1 = vettore delle quote osservate selezionate del bolide (km)
% Velocity = vettore di tutte le velocità osservate del bolide (km/s)
% Velocity1 = vettore delle velocità osservate selezionate del bolide (km/s)
% Tempo1 = tempi osservati selezionati per la posizione (s)
% inclinazione = inclinazione media della traiettoria (radianti)
%
% OUTPUT:
% V_infty = velocità del bolide all'infinito (km/s)
% DV_infty = deviazione standard della velocità all'infinito (km/s)
% fit_velocity = vettore del fit delle velocità con il modello di Ceplecha
% Afin = accelerazione nel punto finale della traiettoria osservata
% (km/s^2)
% A2 = vettore calcolato per l'accelerazione del bolide (km/s^2)
% model_type = modello della velocità utilizzato (Ceplecha o Lineare)
%
% Figure plottate (ma non salvate sull'HDD)
% Plot test della velocità osservate vs quota nel modello di Ceplecha o Lineare
% Plot test delle velocità osservate vs tempo e del fit con il modello di Ceplecha o Lineare
% Plot test dell'accelerazione calcolata del bolide nel modello di Ceplecha o Lineare
%
% Albino Carbognani, versione del 30 aprile 2019


function [V_infty, DV_infty, fit_velocity, Afin, A2, model_type]=Ceplecha_model_velocity2(Height, Height1, Velocity, Velocity1, Tempo1, inclinazione)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stima preliminare della velocità all'infinito del meteoroide con il   %
%   modello empirico di Ceplecha (1961): V(h)=Vinfty+Cv*exp(-k*h)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definizione del modello di Ceplecha per la velocità in funzione della
% quota (Ceplecha, 1961).

Model_Ceplecha_V=@(p, x)(p(1)+p(2)*exp(1).^(-p(3)*x)); 

V2=median(Velocity1); % Valore guess della velocità all'infinito. 

% NOTA: Per il gues della velocità all'infinito si può scegliere fra:
%
% max(Velocity1)
% median(Velocity1)
%
% La funzione mediana è poco sensibile ai punti più scatterati della velocità 
% e a quelli che si trovano verso la fine della traiettoria del bolide
% quando la decelerazione è elevata.

% Sopprime i warning nel caso di fit scadente
% warning('off', 'MATLAB:rankDeficientMatrix')
% warning('off', 'stats:nlinfit:ModelConstantWRTParam') 
% warning('off', 'stats:nlinfit:IterationLimitExceeded')

disp('   ')
disp(' From Ceplecha model velocity2 function: COMPUTE CEPLECHA VELOCITY MODEL  ')
disp('   ')

startingVals = [V2 mean(Velocity1) 0.05];         % Valori iniziali stimati dei parametri p(1), p(2) e p(3)
[coefEsts, R1, J1, CovB1, MSE1] = nlinfit(Height1, Velocity1, Model_Ceplecha_V, startingVals);

% Legge l'eventuale warning che arriva dal bad fit del modello di Ceplecha
[warnmsg, msgid] = lastwarn;                                

% Confronto fra due stringhe. Se sono uguali tf=1, se sono diverse tf=0.
tf1=strcmp(msgid,'MATLAB:rankDeficientMatrix');
tf2=strcmp(msgid,'stats:nlinfit:ModelConstantWRTParam');
tf3 = strcmp(msgid,'stats:nlinfit:IterationLimitExceeded');  
tf4=strcmp(msgid,'stats:nlinfit:IllConditionedJacobian');

% Se il fit con il modello di Ceplecha è buono si prosegue, altrimenti si passa al modello lineare
if (tf1 == 0) && (tf2 == 0) && (tf3 == 0) && (tf4 == 0)

disp('   ')
disp(' From Ceplecha model velocity2 function: CEPLECHA FIT IS GOOD!  ')
disp('   ')    
    
fit_velocity=Model_Ceplecha_V(coefEsts, Height1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stima velocità all'infinito %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_infty=coefEsts(1); 

DV_infty=sqrt(CovB1(1,1)); % Incertezza velocità all'infinito

disp('   ')
disp(' From Ceplecha model velocity2 function: COMPUTE RAW ACCELERATION USING CEPLECHA VELOCITY FIT VS TIME  ')
disp('   ')

% Calcolo del vettore dell'accelerazione con il modello della velocità di Ceplecha

Ccep=coefEsts(2); % Secondo coefficiente modello di Ceplecha
Kcep=coefEsts(3); % Terzo coefficiente modello di Ceplecha

Vh=-(fit_velocity)*sin(inclinazione);                                % Vettore della velocità di caduta verticale (km/s)
A_obs=gradient(Velocity1, Tempo1);                                   % Accelerazione osservata (km/s^2)
A2=-(Ccep*Kcep)*(exp(1).^-(Kcep*Height1)).*Vh;                       % Accelerazione nel modello di Ceplecha in km/s^2
Afin=A2(length(Tempo1));                                             % Accelerazione finale del bolide nel modello di Ceplecha (km/s^2)

model_type='Ceplecha';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figure con il modello delle velocità di Ceplecha (1961) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot quota del bolide vs. velocità 
figure
plot(Height, Velocity, 'r.', 'MarkerSize', 16)          % Valori osservati della velocità
hold on
plot(Height1, fit_velocity, 'k-', 'LineWidth', 3)       % Fit delle velocità 
grid
legend('Observed values','Ceplecha model')
ylabel('Velocity (km/s)','FontSize',20)
xlabel('Height (km)','FontSize',20)
title('Test plot: Ceplecha model - Velocity vs height','FontSize',20)
hold off

% Plot della velocità del bolide vs. tempo 
figure
plot(Tempo1, Velocity1, 'r.', 'MarkerSize', 16)    % Plot delle velocità osservate
hold on
plot(Tempo1, fit_velocity, 'k.', 'MarkerSize', 16) % Plot del fit 
grid
legend('Observed values','Ceplecha model')
xlabel('Time (s)','FontSize',20)
ylabel('Velocity (km/s)','FontSize',20)
title('Test plot: Ceplecha model - Velocity vs time','FontSize',20)
hold off

% Plot dell'accelerazione del bolide vs. tempo
figure
plot(Tempo1, A2, 'k.', 'MarkerSize', 16) 
hold on
plot(Tempo1, A_obs, 'r.', 'MarkerSize', 16)
grid
legend('Ceplecha model', 'Observed values')
xlabel('Time (s)','FontSize',20)
ylabel('Acceleration (km/s^2)','FontSize',20)
title('Test plot: Ceplecha model - Acceleration vs time','FontSize',20)
hold off

% Plot dell'accelerazione del bolide vs. quota
figure
plot(Height1, A2, 'k.', 'MarkerSize', 16) 
hold on
plot(Height1, A_obs, 'r.', 'MarkerSize', 16)
grid
legend('Ceplecha model', 'Observed values')
xlabel('Height (km)','FontSize',20)
ylabel('Acceleration (km/s^2)','FontSize',20)
title('Test plot: Ceplecha model - Acceleration vs height','FontSize',20)
hold off

else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modello lineare fra velocità e quota, adatto per bolidi molto veloci   %
% che si disintegrano subito.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % NOTA: il modello lineare entra in gioco se il fit del modello della velocità di Ceplecha non converge 
   
   clear fit_velocity
   
   disp('   ')
   disp(' From Ceplecha model velocity2 function: BAD CEPLECHA MODEL, SKYP TO LINEAR VELOCITY MODEL  ')
   disp('   ')
    
   Model_Ceplecha_V2=@(p, x)(p(1)+p(2)*x);     % Modello lineare fra velocità e quota x 
   startingVals2 = [V2 0];                     % Valori iniziali stimati dei parametri p(1), p(2), p(3)
   [coefEsts2, R1, J1, CovB2, MSE1] = nlinfit(Height1, Velocity1, Model_Ceplecha_V2, startingVals2);

   fit_velocity=Model_Ceplecha_V2(coefEsts2, Height1);
   
   disp('   ')
   disp(' From Ceplecha model velocity2 function: COMPUTE RAW ACCELERATION USING LINEAR VELOCITY FIT VS TIME  ')
   disp('   ')
   
   V_infty=Model_Ceplecha_V2(coefEsts2, 120);              % La velocità all'infinito del meteoroide è quella a 120 km di quota
   DV_infty=sqrt(CovB2(1,1));                              % Incertezza velocità all'infinito
   Vh=-(fit_velocity)*sin(inclinazione);                   % Vettore della velocità di caduta verticale (km/s)
   A_obs=gradient(Velocity1, Tempo1);                      % Accelerazione osservata (km/s^2)
   A2=coefEsts2(2)*Vh;                                     % Accelerazione nel modello lineare (km/s^2)
   Afin=A2(length(Tempo1));                                % Accelerazione finale del bolide nel modello lineare (km/s^2)
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Plot delle figure con il modello di velocità lineare %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Plot quota del bolide vs. velocità 
   figure
   plot(Height, Velocity, 'r.', 'MarkerSize', 16)          % Valori osservati della velocità
   hold on
   plot(Height1, fit_velocity, 'k-', 'LineWidth', 3)       % Fit della velocità
   grid
   legend('Observed values','Linear model')
   ylabel('Velocity (km/s)','FontSize',20)
   xlabel('Height (km)','FontSize',20)
   title('Test plot: Linear model - Velocity vs height','FontSize',20)
   hold off
   
   % Plot della velocità del bolide vs. tempo 
   figure
   plot(Tempo1, Velocity1, 'r.', 'MarkerSize', 16)    % Plot delle velocità osservate
   hold on
   plot(Tempo1, fit_velocity, 'k.', 'MarkerSize', 16) % Plot del fit 
   grid
   legend('Observed values','Linear model')
   xlabel('Time (s)','FontSize',20)
   ylabel('Velocity (km/s)','FontSize',20)
   title('Test plot: Linear model - Velocity vs time','FontSize',20)
   hold off
   
   % Plot dell'accelerazione del bolide vs. tempo 
   figure
   plot(Tempo1, A2, 'k.', 'MarkerSize', 16)
   hold on
   plot(Tempo1, A_obs, 'r.', 'MarkerSize', 16)
   grid
   legend('Linear model', 'Observed values')
   xlabel('Time (s)','FontSize',20)
   ylabel('Acceleration (km/s^2)','FontSize',20)
   title('Test plot: Linear model - Acceleration vs time','FontSize',20)
   hold off
    
   % Plot dell'accelerazione del bolide vs. quota
   figure
   plot(Height1, A2, 'k.', 'MarkerSize', 16) 
   hold on
   plot(Height1, A_obs, 'r.', 'MarkerSize', 16)
   grid
   legend('Linear model', 'Observed values')
   xlabel('Height (km)','FontSize',20)
   ylabel('Acceleration (km/s^2)','FontSize',20)
   title('Test plot: Linear model - Acceleration vs height','FontSize',20)
   hold off
   
   model_type='  Linear';
   
end

end
