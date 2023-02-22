% MECCANICA CELESTE BARICENTRICA
%
% Funzione per il calcolo degli elementi orbitali baricentrici di un
% meteoroide usando posizione e velocità rispetto al baricentro del 
% Sistema Solare.
%
% BIBLIOGRAFIA:
% L'algoritmo che fa passare dal momento angolare di un corpo celeste agli elementi orbitali è tratto da: 
% T. E. Sterne, An Introduction to Celestial Mechanics, Interscience Publishers, 1960.
%
% Per la posizione e la velocità della Terra rispetto al centro di massa
% del Sistema Solare vedi:
% Hirayama and Kinoshita, 1988. Analytical expressions of the Earth's position and velocity for the calculation of apparent positions. Celestial Mechanics, 41, 389-410.
%
% INPUT:
% jd=giorno giuliano del fireball
% ar, dr=coordinate J2000.0 del radiante vero (corretto per rotazione e
% gravità terrestre)
% v1=velocità geocentrica vera del meteoroide in km/s
% lo=longitudine eclittica del Sole in radianti (J2000.0)
% r_helio=distanza Terra-Sole in UA
%
% OUTPUT:
% lr, br=coordinate eclittiche baricentriche del radiante vero (equinozio
% J2000.0), in radianti.
% ar, dr=coordinate equatoriali baricentriche del radiante vero (equinozio
% J2000.0), in radianti.
% a=semiasse maggiore in UA
% inc=inclinazione dell'orbita in radianti
% go=longitudine del nodo ascendente in radianti
% e=eccentricità
% po=argomento del perielio in radianti
% t=tempo del passaggio al perielio in JD
% v1=velocità baricentrica del meteoroide in UA/d (1UA/d = 1731.456829 km/s)
% v=anomalia vera in radianti
%
% Albino Carbognani (OAVdA)
% Versione del 21 maggio 2019
 
function [lr, br, ar, dr, a, inc, go, e, po, t, v1, v]=meccanica_celeste_baricentrica(jd, ar, dr, v1, lo, r_helio)

% Costanti ausiliarie
C1=57.29577951;                                                             % Conversione da gradi a radianti e viceversa
C2 = 0.000577548;                                                           % conversione km/s <--> UA/d 
ep = 23.4392911/C1;                                                         % obliquita' dell'eclittica al 2000.0 (in gradi ma convertita in radianti)
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

  % calcolo componenti eclittiche della distanza Sole-Terra (UA)

  rx_helio=-r_helio*cos(lo);

  ry_helio=-r_helio*sin(lo);

  % calcolo della posizione eliocentrica del Baricentro del Sistema Solare in AU
  % rispetto all'equinozio medio J2000.0
  
  T=(jd-2451545.0)/36525; % Tempo in secoli giuliani trascorsi dal J2000.0
  
  x_B1=0.004957*sin((124.35 + 3034.91*T)/C1);
  x_B2=0.002718*sin((140.10 + 1222.11*T)/C1);
  x_B3=0.001554*sin((34.36 + 218.49*T)/C1);
  x_B4=0.000835*sin((44.05 + 428.47*T)/C1);
  x_B5=-0.000294;
  x_B6=0.000120*sin((144.4 + 6069.8*T)/C1);
  x_B7=0.000076*sin((95.7 + 2444.2*T)/C1);
  
  y_B1=0.004955*sin((34.36 + 3034.91*T)/C1);
  y_B2=0.002722*sin((50.06 + 1222.11*T)/C1);
  y_B3=0.001554*sin((304.34 + 218.49*T)/C1);
  y_B4=0.000834*sin((314.05 + 428.47*T)/C1);
  y_B5=-0.000339;
  y_B6=0.000120*sin((54.4 + 6069.8*T)/C1);
  y_B7=0.000076*sin((5.7 + 2444.2*T)/C1);
  
  z_B1=0.000118*sin((296.40 + 1222.11*T)/C1);
  z_B2=0.000113*sin((293.89 + 3034.91*T)/C1);
  z_B3=0.000048*sin((172.57 + 218.49*T)/C1);
  
  % posizione eliocentrica del Baricentro
  
  rx_B=x_B1+x_B2+x_B3+x_B4+x_B5+x_B6+x_B7;
  ry_B=y_B1+y_B2+y_B3+y_B4+y_B5+y_B6+y_B7;
  rz_B=z_B1+z_B2+z_B3;
  
  % calcolo componenti eclittiche della distanza Baricentro Sistema
  % Solare-Terra (AU)
  
  rx=rx_helio-rx_B;
  ry=ry_helio-ry_B;
  rz=-rz_B;
  
  r=sqrt(rx^2+ry^2+rz^2); % distanza baricentrica della Terra (AU)
  
  % calcolo componenti eclittiche della velocità baricentrica della Terra
  % (AU/d)
  
  vex_B1=0.01719914*sin((280.46732 + 35999.37285*T)/C1);
  vex_B2=0.00028736*sin((277.996 + 71998.746*T)/C1);
  vex_B3=0.00000715*sin((34.35 + 3034.91*T)/C1);
  vex_B4=0.00000540*sin((275.5 + 107998.1*T)/C1);
  vex_B5=0.00000159*sin((50.1 + 1222.1*T)/C1);
  vex_B6=0.00000039*sin((353.3 + 958465.4*T)/C1);
  vex_B7=0.00000034*sin((52 + 6070*T)/C1);
  vex_B8=0.00000031*sin((168 + 68964*T)/C1);
  vex_B9=0.00000029*sin((172 + 36020*T)/C1);
  vex_B10=0.00000029*sin((209 + 35979*T)/C1);
  vex_B11=0.00000021*sin((264 + 81036*T)/C1);
  vex_B12=0.00000020*sin((2 + 58518*T)/C1);
  vex_B13=0.00000017*sin((314 + 428*T)/C1);
  vex_B14=0.00000016*sin((304 + 218*T)/C1);
  vex_B15=0.00000016*sin((31 + 29930*T)/C1);
  vex_B16=0.00000011*sin((129 + 39034*T)/C1);
  vex_B17=0.00000011*sin((74 + 344*T)/C1);
  vex_B18=0.00000011*sin((273 + 143997*T)/C1);
  vex_B19=0.00000011*sin((255 + 32964*T)/C1);
  vex_B20=0.00000010*sin((53 + 101928*T)/C1);
  vex_B21=0.00000177*T*sin((164 + 71999*T)/C1);
  
  vey_B1=0.01720022*sin((190.46558 + 35999.37285*T)/C1);
  vey_B2=0.00028737*sin((187.995 + 71998.746*T)/C1);
  vey_B3=0.00000715*sin((128.32 + 481266.48*T)/C1);
  vey_B4=0.00000715*sin((304.36 + 3034.91*T)/C1);
  vey_B5=0.00000540*sin((185.5 + 107998.1*T)/C1);
  vey_B6=0.00000159*sin((320.1 + 1222.1*T)/C1);
  vey_B7=0.00000039*sin((263.3 + 958465.4*T)/C1);
  vey_B8=0.00000034*sin((322 + 6070*T)/C1);
  vey_B9=0.00000031*sin((78 + 68964*T)/C1);
  vey_B10=0.00000029*sin((172 + 36020*T)/C1);
  vey_B11=0.00000029*sin((119 + 35979*T)/C1);
  vey_B12=0.00000021*sin((174 + 81036*T)/C1);
  vey_B13=0.00000020*sin((272 + 58518*T)/C1);
  vey_B14=0.00000017*sin((224 + 428*T)/C1);
  vey_B15=0.00000016*sin((214 + 218*T)/C1);
  vey_B16=0.00000016*sin((122 + 29939*T)/C1);
  vey_B17=0.00000011*sin((39 + 39034*T)/C1);
  vey_B18=0.00000011*sin((344 + 344*T)/C1);
  vey_B19=0.00000011*sin((183 + 143997*T)/C1);
  vey_B20=0.00000010*sin((165 + 32964*T)/C1);
  vey_B21=0.00000010*sin((323 + 101928*T)/C1);
  vey_B22=0.00000177*T*sin((74 + 71999*T)/C1);
  
  vez_B1=0.00000065*sin((3.3 + 483202.0*T)/C1);
  vez_B2=0.00000016*sin((204 + 3035*T)/C1);
  vez_B3=0.00000392*T*sin((15.6 + 35999.4*T)/C1);
  
  % componenti eclittiche baricentriche della velocità della Terra, UA/d  
  
  vex=vex_B1+vex_B2+vex_B3+vex_B4+vex_B5+vex_B6+vex_B7+vex_B8+vex_B9+vex_B10+vex_B11+vex_B12+vex_B13+vex_B14+vex_B15+vex_B16+vex_B17+vex_B18+vex_B19+vex_B20+vex_B21;
  vey=vey_B1+vey_B2+vey_B3+vey_B4+vey_B5+vey_B6+vey_B7+vey_B8+vey_B9+vey_B10+vey_B11+vey_B12+vey_B13+vey_B14+vey_B15+vey_B16+vey_B17+vey_B18+vey_B19+vey_B20+vey_B21+vey_B22;
  vez=vez_B1+vez_B2+vez_B3;
  
  ve=sqrt(vex^2+vey^2+vez^2); % Modulo della velocità eliocentrica della Terra, UA/d
  
 % calcolo velocita' eclittica baricentrica del meteoroide, UA/d

  v1x=C2*v1x+vex;

  v1y=C2*v1y+vey;

  v1z=C2*v1z+vez;

  v1=sqrt(v1x*v1x+v1y*v1y+v1z*v1z);

 % Calcolo coordinate eclittiche baricentriche del radiante vero (equinozio
 % J2000.0). 

  br=asin(-(v1z/v1));

  x=-(v1x/(v1*cos(br)));

  y=-(v1y/(v1*cos(br)));

  lr=atan2(y,x);

  if(lr<0)

    lr=lr+2*pi;
    
  end

 % Calcolo coordinate equatoriali baricentriche del radiante vero.
 % Equinozio 2000. 

  dr=asin(cos(ep)*sin(br)+sin(ep)*cos(br)*sin(lr));

  x=(cos(br)*cos(lr))/(cos(dr));

  y=(cos(ep)*cos(br)*sin(lr)-sin(ep)*sin(br))/(cos(dr));

  ar=atan2(y,x);

  if(ar<0)

      ar=ar+2*pi;
  end

 % calcolo semiasse maggiore, periodo e moto medio baricentriche

  x=(2/r)-(v1*v1/U);

  a=1/x;

 % calcolo momento angolare baricentrico del meteoroide per unita' di massa 

  mx=ry*v1z-rz*v1y;

  my=rz*v1x-rx*v1z;

  mz=rx*v1y-ry*v1x;

  m=sqrt(mx*mx+my*my+mz*mz);

 % calcolo inclinazione dell'orbita

  inc=acos(mz/m);

 % calcolo longitudine del nodo ascendente */

 go=atan2(mx/m, -my/m);
 
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