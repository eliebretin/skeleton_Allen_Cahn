function [U_domaine] = surface_bord(X1,X2,X3,epsilon,nums)

N = size(X1);



switch nums;
    case 1
pos1_function = @(theta)   (0.6 + 0.25*cos(theta)).*cos(2*theta)   ;
pos2_function = @(theta)   (0.6 + 0.25*cos(theta)).*sin(2*theta)  ;
pos3_function = @(theta)   (0.25*sin(theta))  ;


Ntheta = 300;
h_theta = 2*pi/Ntheta;
theta = linspace(0,2*pi,Ntheta);

 
distb_boule = 100*ones(N);
r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end

U_domaine = 1/2 - tanh(0.5*distb_boule/(epsilon))/2;

 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% truc 
 case 2,

        
theta0 = 3*pi/24;     
l = 2*0.7*sin(theta0);  
H = 1.7;
l2 = 0.7;

pos1_function = @(theta)  0.7*cos((2*pi-2*theta0)*theta+theta0).*(theta>=0).*(theta<1)  + 0.7*cos(-theta0).*(theta>=1).*(theta<2) ...
+ (0.7*cos(-theta0)- (theta-2)*l2 ).*(theta>=2).*(theta<3) + (0.7*cos(-theta0)- l2 ).*(theta>=3).*(theta<4) ...
+ (0.7*cos(-theta0)- l2 ).*(theta>=4).*(theta<5) + (0.7*cos(-theta0)- l2 ).*(theta>=5).*(theta<6) ...
+ (0.7*cos(-theta0)- l2 +  l2*(theta - 6) ).*(theta>=6).*(theta<7) + (0.7*cos(-theta0)).*(theta>=7);     
pos2_function = @(theta)  0.7*sin((2*pi-2*theta0)*theta+theta0).*(theta>=0).*(theta<1)  + 0.7*sin(-theta0).*(theta>=1).*(theta<2) ...
+ (0.7*sin(-theta0)).*(theta>=2).*(theta<3) + (0.7*sin(-theta0)).*(theta>=3).*(theta<4) ...
+ (0.7*sin(-theta0) + (theta-4)*l ).*(theta>=4).*(theta<5)  + (0.7*sin(-theta0) + l ).*(theta>=5).*(theta<6) ...
+ (0.7*sin(-theta0)+ l ).*(theta>=6).*(theta<7)  + (0.7*sin(-theta0)+ l ).*(theta>=7)  ;     
pos3_function = @(theta)  zeros(size(theta)).*(theta>=0).*(theta<1)  +   H/2*(theta-1).*(theta>=1).*(theta<2) ...
 + H/2*(theta>=2).*(theta<3)+ (H/2 - H*(theta-3)  ).*(theta>=3).*(theta<4) ...
+ (H/2 - H).*(theta>=4).*(theta<5) + (-H/2 + H*(theta-5)).*(theta>=5).*(theta<6) ...
 + (H/2).*(theta>=6).*(theta<7) +  (H/2 - H/2*(theta-7)).*(theta>=7).*(theta<8) +  0.*(theta>8) ;          
NN = 200; 
theta = [linspace(0,1,NN),linspace(1,8,2*NN)];
 Ntheta = length(theta);
 
 distb_boule = 100*ones(N);

r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end

     
U_domaine = 1/2 - tanh(0.5*distb_boule/(epsilon))/2;

    
     
     
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mobius_circle
    case 3,
        
   
 alpha = 1.1;    
pos1_function = @(theta)   alpha*(0.45 + 0.25*cos(theta)).*cos(2*theta)   ;
pos2_function = @(theta)   alpha*(0.45 + 0.25*cos(theta)).*sin(2*theta)  ;
pos3_function = @(theta)   alpha*(0.25*sin(theta))  ;


Ntheta = 300;
h_theta = 2*pi/Ntheta;
theta = linspace(0,2*pi,Ntheta);

 
 distb_boule = 100*ones(N);

r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end



pos1_function = @(theta)   alpha*(-0.2 + (0.55).*cos(theta))   ;
pos2_function = @(theta)   alpha*((0.75).*sin(theta))  ;
pos3_function = @(theta)   alpha*(-(0.25*sin(theta)))  ;


Ntheta = 300;
h_theta = 2*pi/Ntheta;
theta = linspace(0,2*pi,Ntheta);

 r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end




             
U_domaine = 1/2 - tanh(0.5*distb_boule/(epsilon))/2;

    
     
        
        
        

    case 4, %%%%%%%%%%%%cube  
 
%    
 
 L = 0.7; 
%  

 pos2_function = @(theta)  L*((theta>=0).*(theta<1).*(-1) + (theta>=1).*(theta<2).*( -1 + (theta-1)*2)...
   +  (theta>=2).*(theta<3).*(1) + (theta>=3).*(theta<4).*(1 - (theta-3)*2) +  (theta>=4).*(theta<5).*(-1) ...
   + (theta>=5).*(theta<6).*( -1 + (theta-5)*2)  +  (theta>=6).*(theta<7).*(1) +  (theta>=7).*(theta<8).*(1 - (theta-7)*2));

 pos1_function = @(theta)  L*((theta>=0).*(theta<1).*(-1) + (theta>=1).*(theta<2).*(-1)...
   +  (theta>=2).*(theta<3).*(-1+ (theta-2)*2) + (theta>=3).*(theta<4).*(1) +  (theta>=4).*(theta<5).*(1) ...
   + (theta>=5).*(theta<6).*(1)  +  (theta>=6).*(theta<7).*(1  - (theta-6)*2) +  (theta>=7).*(theta<8).*(-1));

pos3_function = @(theta)  L*((theta>=0).*(theta<1).*(-1 + (theta-0)*2) + (theta>=1).*(theta<2).*( 1)...
   +  (theta>=2).*(theta<3).*(1) + (theta>=3).*(theta<4).*(1) +  (theta>=4).*(theta<5).*(1 - 2*(theta-4)) ...
   + (theta>=5).*(theta<6).*( -1)  +  (theta>=6).*(theta<7).*(-1) +  (theta>=7).*(theta<8).*(-1));

  Ntheta = 300;
  theta = linspace(0,8,Ntheta);
  h_theta = 8/Ntheta;
%  
distb_boule = 100*ones(N);

r = 0.04;
for i=1:Ntheta-1,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end
  
  
%   
%   
 pos2_function = @(theta)  L*((theta>=0).*(theta<1).*(-1) + (theta>=1).*(theta<2).*( -1 + (theta-1)*2)...
   +  (theta>=2).*(theta<3).*(1) + (theta>=3).*(theta<4).*(1 - (theta-3)*2) +  (theta>=4).*(theta<5).*(-1) ...
   + (theta>=5).*(theta<6).*( -1 + (theta-5)*2)  +  (theta>=6).*(theta<7).*(1) +  (theta>=7).*(theta<8).*(1 - (theta-7)*2));

 pos3_function = @(theta)  L*((theta>=0).*(theta<1).*(-1) + (theta>=1).*(theta<2).*(-1)...
   +  (theta>=2).*(theta<3).*(-1+ (theta-2)*2) + (theta>=3).*(theta<4).*(1) +  (theta>=4).*(theta<5).*(1) ...
   + (theta>=5).*(theta<6).*(1)  +  (theta>=6).*(theta<7).*(1  - (theta-6)*2) +  (theta>=7).*(theta<8).*(-1));

pos1_function = @(theta)  L*((theta>=0).*(theta<1).*(-1 + (theta-0)*2) + (theta>=1).*(theta<2).*( 1)...
   +  (theta>=2).*(theta<3).*(1) + (theta>=3).*(theta<4).*(1) +  (theta>=4).*(theta<5).*(1 - 2*(theta-4)) ...
   + (theta>=5).*(theta<6).*( -1)  +  (theta>=6).*(theta<7).*(-1) +  (theta>=7).*(theta<8).*(-1));

% 
for i=1:Ntheta-1,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end
  
U_domaine = 1/2 - tanh(0.5*distb_boule/(epsilon))/2;


  
  case 5,
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% pyramide %%%%%%%%%%%%%%%
  theta0 = pi/6;
 p0 = [1,0,-0.5];
  p1 = [-sin(theta0),cos(theta0),-0.5];
  p2 = [-sin(theta0),-cos(theta0),-0.5];
  p3 = [0,0,+0.9];
  
    
 L =  0.9;

 pos1_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(1)*(theta) + p1(1)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p1(1)*(theta-1) +  p3(1)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p3(1)*(theta-2) +  p2(1)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p2(1)*(theta-3) +  p0(1)*(4- theta)  ))  ;

 pos2_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(2)*(theta) + p1(2)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p1(2)*(theta-1) +  p3(2)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p3(2)*(theta-2) +  p2(2)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p2(2)*(theta-3) +  p0(2)*(4- theta)  ))  ;

 pos3_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(3)*(theta) + p1(3)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p1(3)*(theta-1) +  p3(3)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p3(3)*(theta-2) +  p2(3)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p2(3)*(theta-3) +  p0(3)*(4- theta)  ))  ;

h = 10^(-7);
der1_function = @(theta) 1.3*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) 1.3*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) 1.3*(pos3_function(theta + h) -pos3_function(theta))/h  ;


Ntheta = 1000;
theta = linspace(0,4,Ntheta);
  h_theta = 4/Ntheta;
%  
  for n=1:Ntheta, 
  F1 = F1 + exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
  F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
  F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
  end

  
 
 pos1_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(1)*(theta) + p3(1)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p3(1)*(theta-1) +  p2(1)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p2(1)*(theta-2) +  p1(1)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p1(1)*(theta-3) +  p0(1)*(4- theta)  ))  ;

 pos2_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(2)*(theta) + p3(2)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p3(2)*(theta-1) +  p2(2)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p2(2)*(theta-2) +  p1(2)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p1(2)*(theta-3) +  p0(2)*(4- theta)  ))  ;

 pos3_function = @(theta)  L*((theta>=0).*(theta<1).*(p0(3)*(theta) + p3(3)*(1 - theta) ) ...
     + (theta>=1).*(theta<2).*( p3(3)*(theta-1) +  p2(3)*(2- theta) )...
   +  (theta>=2).*(theta<3).*(  p2(3)*(theta-2) +  p1(3)*(3- theta)   ) ...
  + (theta>=3).*(theta<4).*(p1(3)*(theta-3) +  p0(3)*(4- theta)  ))  ;

h = 10^(-7);
der1_function = @(theta) 0.8*1i*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) 0.8*1i*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) 0.8*1i*(pos3_function(theta + h) -pos3_function(theta))/h  ;


Ntheta = 1000;
theta = linspace(0,4,Ntheta);
  h_theta = 4/Ntheta;
%  
  for n=1:Ntheta, 
  F1 = F1 + exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
  F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
  F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
  end

  F1 = F1/1.5;
  F2 = F2/1.5;
  F3 = F3/1.5;
  
  
  %%%%%%%  trois tube 
    
  case 6,
 
  % 
pos1_function = @(theta)  6*(0.4*ones(size(theta))).*(theta.*(1 - theta))/0.9;
pos2_function = @(theta)  6*(-0*0.4*sqrt(3)/2*ones(size(theta))).*(theta.*(1 - theta));
pos3_function = @(theta)  0.9*(-1 + 2*theta) ;

h = 10^(-7);
der1_function = @(theta) (1).*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) (1).*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) (1).*(pos3_function(theta + h) -pos3_function(theta))/h  ;

 Ntheta = 300;
 LL = 1;
 theta = linspace(0,LL,Ntheta);
 h_theta = LL/Ntheta;
 
 distb_boule = 100*ones(N);

r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end

       
  % 
pos1_function = @(theta)  6*(-0.4*1/2*ones(size(theta))).*(theta.*(1 - theta))/0.9;
pos2_function = @(theta)  6*(0.4*sqrt(3)/2*ones(size(theta))).*(theta.*(1 - theta))/0.9;
pos3_function = @(theta)  (1 - 2*theta)*0.9 ;

h = 10^(-7);
der1_function = @(theta) (1/2 + sqrt(3)/2*1i).*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) (1/2 + sqrt(3)/2*1i).*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) (1/2 + sqrt(3)/2*1i).*(pos3_function(theta + h) -pos3_function(theta))/h  ;

 Ntheta = 300;
 LL = 1;
 theta = linspace(0,LL,Ntheta);
 h_theta = LL/Ntheta;
 

r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end

    % 
pos1_function = @(theta)  6*(-0.4*1/2*ones(size(theta))).*(theta.*(1 - theta))/0.9;
pos2_function = @(theta)  6*(-0.4*sqrt(3)/2*ones(size(theta))).*(theta.*(1 - theta))/0.9;
pos3_function = @(theta)  (1 - 2*theta)*0.9 ;

h = 10^(-7);
der1_function = @(theta) (1).*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) (1).*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) (1).*(pos3_function(theta + h) -pos3_function(theta))/h  ;



 Ntheta = 300;
 LL = 1;
 theta = linspace(0,LL,Ntheta);
 h_theta = LL/Ntheta;
 

r = 0.02;
for i=1:Ntheta,
    i
distb_boule =  min(distb_boule,sqrt( (X1- pos1_function(theta(i))/2).^2 + (X2 - pos2_function(theta(i))/2).^2 + (X3  - pos3_function(theta(i))/2).^2) - r);
end
  


     
U_domaine = 1/2 - tanh(0.5*distb_boule/(epsilon))/2;

  
    case 7, %Helice

pos1_function = @(theta)  0.4*cos(4*theta);
pos2_function = @(theta)  0.4*sin(4*theta);
pos3_function = @(theta)  -1 + theta/(2*pi)*2  ;


h = 10^(-7);
der1_function = @(theta) (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
 theta = linspace(0,2*pi,Ntheta);
 h_theta = 2*pi/Ntheta;
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
 

pos1_function = @(theta)  0*cos(6*theta);
pos2_function = @(theta)  0*sin(6*theta);
pos3_function = @(theta)  1 - theta/(2*pi)*2  ;


h = 10^(-7);
der1_function = @(theta) (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 100;
 theta = linspace(0,2*pi,Ntheta);
 h_theta = 2*pi/Ntheta;
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
 
  
  
    case 8,
       
R = 0.75;        
pos1_function = @(theta)   R*cos(theta)   ;
pos2_function = @(theta)   R*sin(theta)  ;
pos3_function = @(theta)   (0.14)*(theta.^0) ;

h = 10^(-7);
der1_function = @(theta)  (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
  h_theta = 2*pi/Ntheta;
 theta = linspace(0,2*pi,Ntheta);

 clf
 plot3(pos1_function(theta),pos2_function(theta),pos3_function(theta))
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
 
  
  pos3_function = @(theta)   -(0.14)*(theta.^0) ;

h = 10^(-7);
der1_function = @(theta)  -(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  -(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  -(pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
  h_theta = 2*pi/Ntheta;
 theta = linspace(0,2*pi,Ntheta);

 clf
 plot3(pos1_function(theta),pos2_function(theta),pos3_function(theta))
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
  
  
  F1 = 1*F1;
  F2 = 1*F2;
  F3 = 1*F3;
  
  case 9,  %%%%%%3 cercle
       
 R = 0.65;        
pos1_function = @(theta)   R*cos(theta)   ;
pos2_function = @(theta)   R*sin(theta)  ;
pos3_function = @(theta)   (0)*(theta.^0) ;

h = 10^(-7);
der1_function = @(theta)  (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
  h_theta = 2*pi/Ntheta;
 theta = linspace(0,2*pi,Ntheta);

 clf
 plot3(pos1_function(theta),pos2_function(theta),pos3_function(theta))
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
        
  R = 0.5;        
pos1_function = @(theta)   (0)*(theta.^0)  ;
pos2_function = @(theta)   R*cos(theta)  ;
pos3_function = @(theta)    R*sin(theta) ;

h = 10^(-7);
der1_function = @(theta)  (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
  h_theta = 2*pi/Ntheta;
 theta = linspace(0,2*pi,Ntheta);

 clf
 plot3(pos1_function(theta),pos2_function(theta),pos3_function(theta))
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
        
 
         
  R = 0.8;        
pos1_function = @(theta)     R*sin(theta) ;
pos2_function = @(theta)     (0)*(theta.^0);
pos3_function = @(theta)    R*cos(theta) ;

h = 10^(-7);
der1_function = @(theta)  (pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta)  (pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta)  (pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
  h_theta = 2*pi/Ntheta;
 theta = linspace(0,2*pi,Ntheta);

 clf
 plot3(pos1_function(theta),pos2_function(theta),pos3_function(theta))
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end
 
 
    otherwise     
        
        
        
        
        
        
% 
pos1_function = @(theta)  0.7*cos(theta);
pos2_function = @(theta)  0.7*sin(theta);
pos3_function = @(theta)  0*0.4*cos(theta/2) ;


h = 10^(-7);
der1_function = @(theta) 1.*(pos1_function(theta + h) -pos1_function(theta))/h  ;
der2_function = @(theta) 1.*(pos2_function(theta + h) -pos2_function(theta))/h  ;
der3_function = @(theta) 1.*(pos3_function(theta + h) -pos3_function(theta))/h  ;


 Ntheta = 300;
 LL = 3*pi/2;
 theta = linspace(0,LL,Ntheta);
 h_theta = LL/Ntheta;
 

 for n=1:Ntheta, 
 F1 =  F1 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der1_function(theta(n))*h_theta  ;
 F2 = F2 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der2_function(theta(n))*h_theta  ; 
 F3 = F3 +  exp(-1*pi*( sqrt(abs(X1-pos1_function(theta(n))).^2 + abs(X2-pos2_function(theta(n))).^2 + abs(X3-pos3_function(theta(n))).^2)/epsilon)).*der3_function(theta(n))*h_theta  ;
 end

 
 
   
  F1 = F1;
  F2 = F2;
  F3 = F3;

end












