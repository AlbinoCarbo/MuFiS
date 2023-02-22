% Funzione per il calcolo delle coordinate del punto di intersezione fra una retta e un piano.
%
% INPUT:
% a, b, c, d = parametri del piano
% XM, YM, ZM = punto da cui passa la retta
% csi, eta, zeta = coseni direttori della retta
%
% OUTPUT:
% x, y, z = coordinate cartesiane del punto di intersezione
%
% Versione del 4 maggio 2018

function [x, y, z]=inter_line_plane(a, b, c, d, XM, YM, ZM, csi, eta, zeta)

 t=(-d-a*XM-b*YM-c*ZM)/(a*csi+b*eta+c*zeta);
 
 x=XM+csi*t;
 y=YM+eta*t;
 z=ZM+zeta*t;
 

