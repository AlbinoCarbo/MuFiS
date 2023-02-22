% MECCANICA CELESTE ELIOCENTRICA
%
% Funzione per il calcolo degli elementi orbitali eliocentrici di un
% meteoroide usando posizione e velocità eliocentriche.
%
% BIBLIOGRAFIA:
% L'algoritmo che fa passare dal momento angolare di un corpo celeste agli elementi orbitali è tratto da: 
% T. E. Sterne, An Introduction to Celestial Mechanics, Interscience Publishers, 1960.
%
% INPUT:
% jd=giorno giuliano del fireball
% ar, dr=coordinate J2000.0 del radiante vero (corretto per rotazione e
% gravità terrestre)
% v1=velocità geocentrica vera del meteoroide in km/s
% lo=longitudine eclittica del Sole in radianti (J2000.0)
% r=distanza Terra-Sole in UA
%
% OUTPUT:
% lr, br=coordinate eclittiche eliocentriche del radiante vero (equinozio
% J2000.0), in radianti.
% ar, dr=coordinate equatoriali eliocentriche del radiante vero (equinozio
% J2000.0), in radianti.
% a=semiasse maggiore in UA
% inc=inclinazione dell'orbita in radianti
% go=longitudine del nodo ascendente in radianti
% e=eccentricità
% po=argomento del perielio in radianti
% t=tempo del passaggio al perielio in JD
% v1=velocità eliocentrica del meteoroide in UA/d (1UA/d = 1731.456829 km/s)
% v=anomalia vera in radianti
%
% Albino Carbognani (OAVdA)
% Versione del 20 maggio 2019
 
function [lr, br, ar, dr, a, inc, go, e, po, t, v1, v]=meccanica_celeste_eliocentrica(jd, ar, dr, v1, lo, r)

% Costanti ausiliarie
C1=57.29577951;                                                             % Conversione da gradi a radianti e viceversa
C2 = 0.000577548;                                                           % conversione km/s <--> UA/d 
ep = 23.4392911/C1;                                                         % obliquita' dell'eclittica al 2000.0 (in gradi ma convertita in radianti)
VE0 = 29.7846918;                                                           % velocita' orbitale media della Terra, km/s
U = 0.000295912;                                                            % costante di gravitazione
 
  % Calcolo coordinate eclittiche geocentriche del radiante (corrette per la rotazione
  % e la gravita' terrestre ma non per la velocita' orbitale). Equinozio
  % J2000.0.

  br=asin(cos(ep)*sin(dr)-sin(ep)*cos(dr)*sin(ar));

  x=cos(dr)*cos(ar)/cos(br);

  y=(sin(dr)*sin(ep)+cos(dr)*cos(ep)*sin(ar))/cos(br);

  lr=atan2(y,x);

  if(lr<0)

    lr=lr+2*pi;
    
  end
  
  % calcolo velocità eclittica geocentrica del meteoroide

  v1x=-v1*cos(br)*cos(lr);

  v1y=-v1*cos(br)*sin(lr);

  v1z=-v1*sin(br);

  % calcolo componenti eclittiche della distanza Sole-Terra

  rx=-r*cos(lo);

  ry=-r*sin(lo);

 % calcolo componenti eclittiche della velocità della Terra   
  
  ve=VE0*sqrt((2/r)-1); % Modulo della velocità eliocentrica della Terra

  vex=ve*sin(lo+0.01672*sin(lo-1.796));

  vey=-ve*cos(lo+0.01672*sin(lo-1.796));

 % calcolo velocita' eclittica eliocentrica del meteoroide, UA/d

  v1x=C2*(v1x+vex);

  v1y=C2*(v1y+vey);

  v1z=C2*v1z;

  v1=sqrt(v1x*v1x+v1y*v1y+v1z*v1z);

 % Calcolo coordinate eclittiche eliocentriche del radiante vero (equinozio
 % J2000.0). 

  br=asin(-(v1z/v1));

  x=-(v1x/(v1*cos(br)));

  y=-(v1y/(v1*cos(br)));

  lr=atan2(y,x);

  if(lr<0)

    lr=lr+2*pi;
    
  end

 % Calcolo coordinate equatoriali eliocentriche del radiante vero.
 % Equinozio 2000. 

  dr=asin(cos(ep)*sin(br)+sin(ep)*cos(br)*sin(lr));

  x=(cos(br)*cos(lr))/(cos(dr));

  y=(cos(ep)*cos(br)*sin(lr)-sin(ep)*sin(br))/(cos(dr));

  ar=atan2(y,x);

  if(ar<0)

      ar=ar+2*pi;
  end

 % calcolo semiasse maggiore, periodo e moto medio

  x=(2/r)-(v1*v1/U);

  a=1/x;

 % calcolo momento angolare del meteoroide per unita' di massa 

  mx=ry*v1z;

  my=-rx*v1z;

  mz=rx*v1y-ry*v1x;

  m=sqrt(mx*mx+my*my+mz*mz);

 % calcolo inclinazione dell'orbita

  inc=acos(mz/m);

 % calcolo longitudine del nodo ascendente */


  if(br>0)

    go=lo;
    
  end

  if(br<0)

    go=lo-pi;
    
  end

  if(go<0)

     go=go+2*pi;
     
  end

 % calcolo eccentricità

  e=sqrt(1-(m*m)/(U*a));

 % calcolo anomalia vera

  vr=(v1x*rx+v1y*ry)/r;

  x=((m*m)/(r*U)-1)/e;

  y=(m*vr)/(U*e);

  v=atan2(y,x);

  if(v<0)

    v=v+2*pi;
    
  end

 % calcolo argomento del perielio

  if(br>0)

    po=pi-v;
    
  end

  if(br<0)

    po=-v;
    
  end

  if(po<0)

    po=po+2*pi;
    
  end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % calcoli per l'orbita ellittica  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(e<1)

     % calcolo periodo e moto medio

     n=sqrt(U/(a*a*a));

     p=2*pi/n;

     % calcolo anomalia eccentrica 

     x=(1-(r/a))/e;

     y=(sin(v)*(1-e*x))/(sqrt(1-e*e));

     ae=atan2(y,x);

     if(ae<0)

        ae=ae+2*pi;
        
     end

    % calcolo anomalia media

     am=ae-e*sin(ae);

    if(am<0)

       am=am+2*pi;
       
    end

   % calcolo data del passaggio al perielio

    t=jd-(am/n);  

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calcolo del passaggio al perielio per l'orbita parabolica %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(e==1)
      
    q1=(m*m)/(2*U); z1=tan(v/2);

    t=jd-sqrt(2*(q1*q1*q1)/U)*(z1-(z1*z1*z1)/3);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calcolo del passaggio al perielio per l'orbita iperbolica %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(e>1)

    f=sqrt((e-1)/(e+1))*tan(v/2);

    f=2*log(sqrt((f+1)/(1-f)));

    t=jd-sqrt(-(a*a*a)/U)*(e*sinh(f)-f);

  end