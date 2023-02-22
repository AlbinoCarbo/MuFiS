% Fireball progenitor/Monte Carlo clones and Solar System orbits plotter
% 
% Calcola, plotta e salva una figura con le orbite - in scala e con la corretta posizione - di tutti i pianeti 
% visti dal polo nord dell'Eclittica. La posizione dei pianeti è quella al momento della caduta del bolide.
% Viene anche tracciata, tratteggiata, l'orbita del pianeta nano Plutone.
%
% Alle orbite ellittiche planetarie viene sovrapposta l'orbita nominale (ellittica o iperbolica) 
% del meteoroide progenitore del bolide più quelle (eventuali, dipende dall'input) dei cloni Monte Carlo 
% supponendo l'inclinazione orbitale nulla. 
% 
% Attenzione che non tutte le orbite planetarie potrebbero essere visibili, dipende dall'orbita del progenitore del bolide.
% Infatti, la figura viene automaticamente ridimensionata in modo che contenga tutta l'orbita nominale del fireball
% e quelle degli eventuali cloni. Questa funzione fa anche un plot a punti random della Fascia Principale degli asteroidi 
% (fra 2.2 e 3.6 UA) e della Fascia di Kuiper (fra 30 e 55 UA).
%
% NOTA: Per visualizzare il Sistema Solare fino alla Fascia di Kuiper basta cambiare i limiti degli assi 
% direttamente dai comandi della figura visualizzata sullo schermo.
%
% INPUT:
% fireball_name = nome del fireball (stringa)
% working_path = path dove salvare l'immagine delle orbite
% jd = giorno giuliano del fireball
% a = numero oppure array con i semiassi maggiori delle orbite del meteoroide progenitore più quelli dei cloni Monte Carlo (UA)
% e = numero oppure array con l'eccentricità delle orbite del meteoroide progenitore più quelle dei cloni Monte Carlo
% omega = numero oppure array con la longitudine del perielio J2000.0 delle orbite.
% del meteoroide progenitore più quelle dei cloni Monte Carlo. 
%
% NOTA: longitudine del perielio = longitudine nodo ascendente + argomento del perielio (gradi)
% N = numero di cloni Monte Carlo
%
% OUTPUT:
% Figura con il plot delle orbite dei pianeti/pianeti nani e del meteoroide
% progenitore.
%
% BIBLIOGRAFIA: per il calcolo delle posizioni dei pianeti sono stati usati
% gli algosritmi del libro di J. Meeus, "Astronomical Algorithms",
% Willmann-Bell, 1991 (prima edizione), capitoli 29 e 30.
%
% Per la formula del seno e coseno dell'anomalia vera è stato consultato il
% libro di Craig A. Kluever, Space Flight Dynamics, Wiley & Sons, 2018, formule 4.13 e 4.29
%
% Albino Carbognani, versione del 1 marzo 2019

function []=Orbit_plotter_monte_carlo(fireball_name, working_path, jd, a, e, omega, N)

% Definizione delle costanti e delle funzioni

C=57.29577951;              % Conversione da gradi a radianti e viceversa
T=(jd-2451545)/36525;       % Tempo in secoli giuliani a partire dal 1.5 gennaio 2000
eps=0.000017455;            % Precisione per la soluzione dell'equazione di Keplero (radianti)
v_orbita=0:0.1:360;         % Valori del vettore dell'anomalia vera per cui disegnare l'orbita ellittica (gradi)

% Nome del file/titolo dell'immagine
if N == 1
  string=[fireball_name ' - Solar System and nominal orbit diagram '];      % Per la figura con orbita nominale
else
  string=[fireball_name ' - Solar System and Monte Carlo orbits diagram ']; % Per la figura con i cloni Monte Carlo  
end

% Defiizione delle funzioni polari per raggio vettore, seno e coseno dell'anomalia
% vera

raggio_vettore=@(x, y, z)((x*(1-y^2))./(1+y*cos(z)));       % Definizione dell'equazione dell'orbita ellittica o iperbolica (x=semiasse maggiore, y=eccentricità, z=anomalia vera). 
                                                            % NOTA: se orbita ellittica x>0 e y<1, se iperbolica x<0 e y>1. In entrambi i casi x(1-y^2) > 0
coseno_v=@(x, y)((cos(x)-y)/(1-y*cos(x)));                  % Coseno anomalia vera: x anomalia eccentrica in radianti, y eccentricità
seno_v=@(x, y)(((sqrt(1-y^2))*sin(x))/(1-y*cos(x)));        % Seno anomalia vera: x anomalia eccentrica in radianti, y eccentricità


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elementi orbitali J2000.0, orbite e posizioni per i pianeti %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=zeros(1, 8);                     % Pre-allocazione del vettore delle anomalie vere per le posizioni dei pianeti
r=zeros(1, 8);                     % Pre-allocazione del vettore delle distanze eliocentriche per le posizioni dei pianeti
longitudine_perielio=zeros(1, 8);

% Mercurio
a1=0.38709893;                                   % Semiasse maggiore (UA)
e1=0.20563069;                                   % Eccentricità
omega1=77.456119;                                % Longitude of perihelion (gradi)
L1=252.250906+149472.6746358*T-0.00000535*T^2;   % Longitudine eclittica J2000.0 (gradi)
M1=L1-omega1;                                    % Anomalia media (gradi)
R1=raggio_vettore(a1, e1, (v_orbita-omega1)/C);  % Orbita completa del pianeta

% Calcolo posizione di Mercurio
E1 = keplerEq(M1/C,e1,eps);                        % Anomalia eccentrica di Mercurio al momento della caduta del bolide (radianti)
v(1)=C*atan2(seno_v(E1, e1), coseno_v(E1, e1));    % Anomalia vera di Mercurio al momento della caduta del bolide (gradi)
r(1)=raggio_vettore(a1, e1, v(1)/C);               % Raggio vettore vero di Mercurio al momento della caduta del bolide (UA)
longitudine_perielio(1)=omega1;

% Venere
a2=0.72333199; 
e2=0.00677323;
omega2=131.563707;
L2=181.979801+58517.8156760*T+0.00000165*T^2;
M2=L2-omega2;
R2=raggio_vettore(a2, e2, (v_orbita-omega2)/C);

% Calcolo posizione di Venere
E2 = keplerEq(M2/C,e2,eps);                        
v(2)=C*atan2(seno_v(E2, e2), coseno_v(E2, e2));    
r(2)=raggio_vettore(a2, e2, v(2)/C);      
longitudine_perielio(2)=omega2;

% Terra
a3=1.00000011;
e3=0.01671022;
omega3=102.937348;
L3=100.466449+35999.3728519*T-0.00000568*T^2;
M3=L3-omega3;
R3=raggio_vettore(a3, e3, (v_orbita-omega3)/C);

% Calcolo posizione della Terra
E3 = keplerEq(M3/C,e3,eps);                        
v(3)=C*atan2(seno_v(E3, e3), coseno_v(E3, e3));    
r(3)=raggio_vettore(a3, e3, v(3)/C);     
longitudine_perielio(3)=omega3;

% Marte
a4=1.52366231;
e4=0.09341233;
omega4=336.060234;
L4=355.433275+19140.2993313*T+0.00000261*T^2;
M4=L4-omega4;
R4=raggio_vettore(a4, e4, (v_orbita-omega4)/C);

% Calcolo posizione di Marte
E4 = keplerEq(M4/C,e4,eps);                        
v(4)=C*atan2(seno_v(E4, e4), coseno_v(E4, e4));    
r(4)=raggio_vettore(a4, e4, v(4)/C); 
longitudine_perielio(4)=omega4;

% Giove
a5=5.20336301;
e5=0.04839266;
omega5=14.331309;
L5=34.351484+3034.9056746*T-0.0008501*T^2;
M5=L5-omega5;
R5=raggio_vettore(a5, e5, (v_orbita-omega5)/C);

% Calcolo posizione di Giove
E5 = keplerEq(M5/C,e5,eps);                        
v(5)=C*atan2(seno_v(E5, e5), coseno_v(E5, e5));    
r(5)=raggio_vettore(a5, e5, v(5)/C); 
longitudine_perielio(5)=omega5;

% Saturno
a6=9.53707032;
e6=0.05415060;
omega6=93.056787;
L6=50.077471+1222.1137943*T+0.00021004*T^2;
M6=L6-omega6;
R6=raggio_vettore(a6, e6, (v_orbita-omega6)/C);

% Calcolo posizione di Saturno
E6 = keplerEq(M6/C,e6,eps);                        
v(6)=C*atan2(seno_v(E6, e6), coseno_v(E6, e6));    
r(6)=raggio_vettore(a6, e6, v(6)/C);  
longitudine_perielio(6)=omega6;

% Urano
a7=19.19126393;
e7=0.04716771;
omega7=173.005159;
L7=314.055005+428.4669983*T-0.00000486*T^2;
M7=L7-omega7;
R7=raggio_vettore(a7, e7, (v_orbita-omega7)/C);

% Calcolo posizione di Urano
E7 = keplerEq(M7/C,e7,eps);                        
v(7)=C*atan2(seno_v(E7, e7), coseno_v(E7, e7));    
r(7)=raggio_vettore(a7, e7, v(7)/C);  
longitudine_perielio(7)=omega7;

% Nettuno
a8=30.06896348;
e8=0.00858587;
omega8=48.123691;
L8=304.348665+218.4862002*T;
M8=L8-omega8;
R8=raggio_vettore(a8, e8, (v_orbita-omega8)/C);

% Calcolo posizione di Nettuno
E8 = keplerEq(M8/C,e8,eps);                        
v(8)=C*atan2(seno_v(E8, e8), coseno_v(E8, e8));    
r(8)=raggio_vettore(a8, e8, v(8)/C);  
longitudine_perielio(8)=omega8;

% Plutone
a9=39.48168677;
e9=0.24880766;
omega9=224.06676;
R9=raggio_vettore(a9, e9, (v_orbita-omega9)/C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dei raggi vettori delle orbite (ellittiche o iperboliche) dei meteoroidi cloni Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vettori dell'anomalia vera (in gradi) per cui disegnare le orbite,
% ellittiche o iperboliche, del meteoroide e i suoi cloni.
%
% NOTA: l'iperbole non viene plottata su tutti i 360° perché altrimenti Matlab
% disegna anche gli asintoti e il secondo ramo dell'iperbole che non ha il fuoco nel Sole.

for i=1:N
    if a(i)>0
        f_orbita{i}=v_orbita;                                                         % Vettore anomalia vera per orbite ellittiche
        
    else
         Theta_infty=C*acos(-1/e(i));                                                 % Angolo fra l'asintoto dell'iperbole e l'asse centrale dell'iperbole (gradi)
         f_orbita{i}=(omega(i)-(Theta_infty-1):0.1:omega(i)+(Theta_infty-1));         % Vettore anomalia vera per le orbite iperboliche

    end
end
  
% Calcolo dei raggi vettori delle orbite ellittiche o iperboliche
% per il meteoroide e i suoi cloni
    
for i=1:N  
      R{i}=raggio_vettore(a(i), e(i), (f_orbita{i}-omega(i))/C);                      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot delle orbite eliocentriche dei pianeti e del meteoroide %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

% Dimensionamento automatico della figura in UA (Sole al centro)
if N == 1
dimension_solar_system=abs((a(1)*(1+e(1))+1)); 

if dimension_solar_system > 50
    dimension_solar_system=50; % Limite superiore alle dimensioni del Sistema Solare plottato
end

xlim([-dimension_solar_system dimension_solar_system]) % Set limiti asse x per contenere l'orbita del bolide
ylim([-dimension_solar_system dimension_solar_system]) % Set limiti asse y per contenere l'orbita del bolide
else
    
dimension_solar_system=max(abs(a.*(1+e)+1));

if dimension_solar_system > 50
    dimension_solar_system=50;  % Limite superiore alle dimensioni del Sistema Solare plottato
end

xlim([-dimension_solar_system dimension_solar_system]) % Set limiti asse x per contenere le orbite dei cloni Monte Carlo
ylim([-dimension_solar_system dimension_solar_system]) % Set limiti asse y per contenere le orbite dei cloni Monte Carlo
end

title(string,'FontSize',17)
ylabel('Distance(AU)','FontSize',17)
xlabel('Distance (AU)','FontSize',17)
daspect([1 1 1]) % Stessa unità di misura sugli assi in modo da non deformare il plot
grid

% Plot della direzione del punto gamma
x0 = [0.84 0.87]; % In unità della dimensione X della figura
y0 = [0.5 0.5];  % In unità della dimensione Y della figura
annotation('textarrow',x0, y0,'String',' \gamma ')

% Plot random degli asteroidi della fascia asteroidale

n = 1000; % Numero di asteroidi da simulare
theta = rand(1,n)*(2*pi); % Anomalia vera asteroidi
r0 = 2.2+1.4*(rand(1,n));  % Raggio vettore asteroidi
x = r0.*cos(theta);
y = r0.*sin(theta);
hold on
plot(x,y,'k.');
hold on

% Plot random degli asteroidi della fascia di Kuiper
n1 = 2000; % Numero di TNO da simulare
theta1 = rand(1,n1)*(2*pi); % Anomalia vera asteroidi
r1 = 30.0+25*(rand(1,n1));  % Raggio vettore asteroidi
x = r1.*cos(theta1);
y = r1.*sin(theta1);
hold on
plot(x,y,'k.');
hold on

% Disegno del Sole

viscircles([0 0], 0.1, 'Color','k'); 
hold on
plot(0, 0, 'y.', 'MarkerSize', 20);
hold on

% Plot dell'orbita, ellittica o iperbolica, del meteoroide progenitore del bolide.
% Se N = 1 si fa il plot dell'orbita nominale del meteoroide progenitore. 
% Se N > 1 si fa il plot anche dei cloni Monte Carlo.

for i=1:N
    
       bolide{i}=polar(f_orbita{i}/C, R{i});         % Orbita ellittica o iperbolica
       set(bolide{i},'LineWidth',2, 'Color','m')
       hold on
end

% Plot dell'orbita nominale del meteoroide al di sopra delle orbite dei cloni Monte Carlo
% con colore diverso in modo da farla risaltare sulle altre.

if N > 1
       
        bolide{1}=polar(f_orbita{1}/C, R{1});         % Orbita ellittica o iperbolica
        set(bolide{1},'LineWidth',2, 'Color','k')
        hold on
end

% Plot orbite pianeti e Plutone

plutone=polar(v_orbita/C, R9);
set(plutone,'LineWidth',2, 'Color','g', 'LineStyle', '--')
hold on
nettuno=polar(v_orbita/C, R8);
set(nettuno,'LineWidth',2, 'Color','b')
hold on
urano=polar(v_orbita/C, R7);
set(urano,'LineWidth',2, 'Color','b')
hold on
saturno=polar(v_orbita/C, R6);
set(saturno,'LineWidth',2, 'Color','k')
hold on
giove=polar(v_orbita/C, R5);
set(giove,'LineWidth',2, 'Color','y')
hold on
marte=polar(v_orbita/C, R4);
set(marte,'LineWidth',2, 'Color','r')
hold on
terra=polar(v_orbita/C, R3);
set(terra,'LineWidth',2, 'Color','b')
hold on
venere=polar(v_orbita/C, R2);
set(venere,'LineWidth',2, 'Color','g')
hold on
mercurio=polar(v_orbita/C, R1);
set(mercurio,'LineWidth',2, 'Color','c')
hold on

% Plot della posizione dei pianeti sull'orbita

for k = 1 : length(v)
    % Trasformazione da coordinate polari a cartesiane con rotazione
    % di un angolo pari alla longitudine del perielio, lo stesso angolo di cui sono state ruotate le orbite.
    [x, y] = pol2cart((v+longitudine_perielio)/C, r); 
    plot(x, y, 'k.', 'MarkerSize', 20);
    hold on
end

hold off

cd(working_path)
saveas(gcf, string, 'bmp')   % Salva la figura in bmp

end