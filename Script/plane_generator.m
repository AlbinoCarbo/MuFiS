% Funzione per il calcolo delle costanti di un piano a partire da due vettori giacenti sul piano e un punto
% da cui il piano deve passare.
%
% INPUT:
% A1, A2, A3 = componenti del primo vettore
% B1, B2, B3 = componenti del secondo vettore
% X, Y, Z = punto da cui deve passare il piano
%
% OUTPUT:
% a, b, c, d = componenti del piano contenente i due vettori e passante per il punto X, Y, Z.
%
% Versione del 5 maggio 2018

function [a, b, c, d]=plane_generator(A1, A2, A3, B1, B2, B3, X, Y, Z)

A=[A1 A2 A3];

B=[B1 B2 B3];

C=cross(A, B)/(norm(cross(A, B)));

a=C(1); b=C(2); c=C(3); 

d=-X*a-Y*b-Z*c; % Condizione del passaggio per il punto X, Y, Z.
