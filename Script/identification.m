% IDENTIFICATION
%
% Parametri formali di input:
% alfa=AR geocentrico del radiante vero in radianti
% delta=Dec geocentrica del radiante vero in radianti
% vg=velocità geocentrica vera del meteoroide in km/s
% epsilon=obliquità dell'Eclittica in radianti
% lambda=longitudine eclittica della Terra in radianti
%
% Output:
% U=modulo velocità geocentrica del meteoroidee espressa in unità della velocità
% media della Terra
% Ux, Uy, Uz=componenti velocità geocentrica secondo questo sistema di riferimento:
% asse z ortogonale all'Eclittica, asse y nella direzione del moto della
% Terra, asse x in direzione opposta al Sole.
%
% theta=angolo in radianti fra il vettore della velocità geocentrica U e la
% direzione y del moto della Terra
% fi=angolo in radianti fra il piano y-z e il piano individuato da U e l'asse y.
% a=semiasse orbitale (UA)
% e=eccentricità orbita
% i=inclinazione orbita (radianti)
%
% Albino Carbognani (OAS)
% Versione del 20 marzo 2020

function [U, Ux, Uy, Uz, theta, fi, aa, ee, ii]=identification(alfa, delta, vg, epsilon, lambda)

Vel_antiradiante=(vg/29.785)*[-cos(delta)*cos(alfa), -cos(delta)*sin(alfa), -sin(delta)];

% Matrice di rotazione attorno all'asse x
RX=[1, 0, 0; 0, cos(-epsilon), -sin(-epsilon); 0, sin(-epsilon), cos(-epsilon)];

% Matrice di rotazione attorno all'asse z
RZ=[cos(-lambda), -sin(-lambda), 0; sin(-lambda), cos(-lambda), 0; 0, 0, 1];

UV=RZ*(RX*Vel_antiradiante');

Ux=UV(1); Uy=UV(2); Uz=UV(3);

U=sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

theta=acos(Uy/U);

fi=atan2(Ux, Uz);

if fi < 0
    fi=fi+2*pi;   % Riduzione di fi ad un valore sempre positivo compreso fra 0 e 2*pi
end

% Calcolo elementi orbitali a, e, i
aa=1/(1-U^2-2*Uy);
ee=sqrt(U^4+4*Uy^2+(Ux^2)*(1-U^2-2*Uy)+4*Uy*(U^2));
ii=atan2(Uz, 1+Uy);
    
end



