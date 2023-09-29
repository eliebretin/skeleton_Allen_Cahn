function U = condition_initiale3d(N,epsilon,no_cond);
U = zeros(N,N,N);
R0 = 0.12;
X1 = zeros(N,N,N);
x= linspace(-1/2,1/2,N);
x1 = ((x)'*ones(1,N));

for i=1:N
X1(:,:,i)=x1;%+ rand(N,N)/(2*N);
end
X2 = permute(X1,[2 1 3]);
X3 = permute(X1,[3 2 1])-0*0.2;



%%%%%%%%%%%%%%%%%%% sphere %%%%%%%%%%%%%%%%%%%
if (no_cond==1),
X1 = 0.9*(X1+0*0.02).^2;
X2 = 0.9*(X2 +0*0.03).^2;
X3 = 0.9*(X3+0.05).^2;
R0 = 0.45;
U=sqrt(X1+X2+X3)-R0;
U = 1/2 - 1/2*tanh(U/(2*epsilon));
%%%%%%%%%%%%%%%%%% sphere troue %%%%%%%%%%%%%%%%%%%%%%
elseif (no_cond==6),
X1_2 = (X1).^2;
X2_2 = X2.^2;
X3_2 = (X3).^2;
R0 = 0.4;
U=X1_2+X2_2+X3_2-R0^2;
U = 1/2 - 1/2*tanh(U/epsilon);

X1 = (X1).^2;
X2 = X2.^2;
X3 = (X3-0.10).^2;
R0 = 0.32;
U_temp=X1+X2+X3-R0^2;
U_temp = 1/2 - 1/2*tanh(U_temp/epsilon/2);
U = min(max(U - U_temp,0),1); 


%%%%%%%%%%%%%%%%%%% tore  %%%%%%%%%%%%%%%%%%%%%
elseif (no_cond==2)
X1 = (X1+0*0.00).^2;
X2 = X2.^2;
X3 = X3.^2;    
R1 = 0.3;
%R0 = 0.15;
R0 = 0.1;
U = sqrt((sqrt(X1 +X2 )-R1).^2 + X3) - 2*R0;
U = 1/2 - 1/2*tanh(U/epsilon/2);


elseif (no_cond==3)
%%%%%%%%%%%%%%%%%%% condition dordel  
U = 0.25+ rand(N,N,N)/75  - X3; 
%U = 0.25  + ((cos(8*pi*X1) + cos(8*pi*X2))/12 + (cos(2*pi*X1) + cos(2*pi*X2))/24)  - X3;
U = 1/2 + sign(U)/2;

%U(1:N/64,:,:) = 0*(0.25 - X3(1:N/64,:,:));
%U(N-N/64:N,:,:) = 0*(0.25  - X3(N-N/64:N,:,:));
%U(:,1:N/64,:) = 0*(0.25 - X3(:,1:N/64,:));
%U(:,N-N/64:N,:) = 0*(0.25  - X3(:,N-N/64:N,:));


U(:,:,1:N/2) = U(:,:,N:-1:N/2+1); 

%%%% carre %%%%%%%%%%%%
elseif (no_cond==4)

  dist = max(max(max(max(max(X1 - 0.4,-0.45 - X1),X2 - 0.45),-0.4 - X2),X3 - 0.45),-0.45 - X3);
  U = 1/2 - 1/2*tanh(dist/epsilon/2);

 
 elseif (no_cond==24)
 dist = max(max(max(max(max(max(X1 - 0.4,-0.45 - X1),X2 - 0.45),-0.4 - X2),X3 - 0.45),-0.45 - X3),0.1-sqrt(X1.^2 + X2.^2));
 
  U = 1/2 - 1/2*tanh(dist/epsilon/2);

  elseif (no_cond==26)
 dist = max(max(max(max(max(max(X1 - 0.4,-0.45 - X1),X2 - 0.45),-0.4 - X2),X3 - 0.45),-0.45 - X3),0.1-sqrt((X2-0.2).^2 + (X3-0.2).^2));
 
  U = 1/2 - 1/2*tanh(dist/epsilon/2);

 
 
% %%%%%%%%%%%%%%%%%%% double torre
elseif (no_cond==5)
X1_1 = (X1+0.125).^2;
X2_1 = X2.^2;
X3_1 = X3.^2;    
R1 = 0.25;
R0 = 0.07;
U = (sqrt(X1_1 +X2_1 )-R1).^2 + X3_1 - R0^2;
U = 1/2 - 1/2*tanh(U/epsilon*10);

X1_2 = (X1-0.125).^2;
X2_2 = X2.^2;
X3_2 = X3.^2;    
R1 = 0.25;
R0 = 0.07;
U_temp = (sqrt(X1_2 +X3_2 )-R1).^2 + X2_2 - R0^2;
U_temp = 1/2 - 1/2*tanh(U_temp/epsilon*10);
U = min(max(U + U_temp,0),1); 
   
%%%%%%%%%%%%%%%%%%%%%% double cube
elseif (no_cond==7)
%U(N/8+1:N/2+N/8,N/8+1:N/2+N/8,N/8+1:N/2+N/8)=ones(N/2,N/2,N/2);
%U(N/2+1-N/8:N/2+N/8,N/2+1-N/8:N/2+N/8,N/2+N/8+1:N/2+N/8+11*N/32)=ones(2*N/8,2*N/8,11*N/32);
%U(:,:,:)=U(N:-1:1,N:-1:1,1:N);

U(N/4-3*N/16+1:N/4+3*N/16,N/2-N/4+1:N/2+N/4,N/2-N/4+1:N/2+N/4)=ones(6*N/16,N/2,N/2);
U(3*N/4-3*N/16+1:3*N/4+3*N/16,N/2-N/4+1:N/2+N/4,N/2-N/4+1:N/2+N/4)=ones(6*N/16,N/2,N/2);

    
%%%%%%%%%%%%%%%%%%%%%% cylindre
elseif (no_cond==9)
R_1 = 0.1;
R_2 = 1;
U_temp = max((sqrt(X1.^2 +X2.^2 )-R_1), X3.^2 - R_2^2);
U = 1/2 - 1/2*tanh(U_temp/epsilon/2);

%%%%%%%%%%%%%%%%%%%%%% cylindre bord
elseif (no_cond==10)
R_1 = 0.4;
R1_l = 0.1;
R_2 = 0.1;
R2_l = 0.05,
U_temp = max(max((sqrt(X1.^2 +X3.^2 )-R_1),(X2.^2 - R_2^2)),-(X2.^2 - (R_2-R2_l)^2));
%U_temp = max(max(( abs(X1) +abs(X3) -R_1),(X2.^2 - R_2^2)),-(X2.^2 - R_3^2));
U = 1/2 - 1/2*sign(U_temp/epsilon/2);
U_temp = max(max((sqrt(X1.^2 +X3.^2 )-(R_1-R1_l)),(X2.^2 - R_2^2)),-(X2.^2 - (R_2-R2_l)^2));
%U_temp = max(max((abs(X1) +abs(X3) -R_1),(X2.^2 - R_2^2)),-(X2.^2 - R_3^2));
U = U - (1/2 - 1/2*sign(U_temp/epsilon/2));


%%%%%%%%%%%%%% double cylindre %%%%%%%%%%%%%%

elseif (no_cond==12)
R_1 = 0.1;
R_2 = 1;
U_temp1 = max((sqrt((X1-0.15).^2 +X3.^2 )-R_1), X2.^2 - R_2^2);

%U_temp2 = max((sqrt((X1+0.15).^2 +X3.^2 )-R_1), X2.^2 - R_2^2);
U_temp2 = max((sqrt((X1+0.15).^2 +X2.^2 )-R_1), X3.^2 - R_2^2);

U = (1/2 - 1/2*tanh(0.5*U_temp1/epsilon)) + (1/2 - 1/2*tanh(0.5*U_temp2/epsilon));

%%%%%%%%%%%%%%% double sphere %%%%%%%%%%%%%%%%%%

elseif (no_cond==11)
R0 = 0.16;
temp1=(X1-0.22).^2+X2.^2+X3.^2-R0^2;
temp2=(X1+0.22).^2+X2.^2+X3.^2-R0^2;
U = 1/2 - 1/2*tanh(temp1/epsilon) + 1/2 - 1/2*tanh(temp2/epsilon)  ;
%%%%%%%%%%%%%%% double sphere %%%%%%%%%%%%%%%%%%
elseif (no_cond==13)
U(1:N/2,:,:) = ones(N/2,N,N);
U(N/2-N/4:N/2-1,N/2-N/4+1:N/2+N/4,N/2-N/4+1:N/2+N/4) = 0*ones(N/4,N/2,N/2);
U(N/2+1:N/2+N/4,N/2-N/4+1:N/2+N/4,N/2-N/4+1:N/2+N/4) = ones(N/4,N/2,N/2);
U(N/2,:,:) = 1/2*ones(1,N,N);

%%%%%%%%%%%%%%% double sphere %%%%%%%%%%%%%%%%%%
elseif (no_cond==14)
U(1:N/2,1:N/2,1:N/2) = ones(N/2,N/2,N/2);
U(1:N/2,N/2+1:N,N/2+1:N) = ones(N/2,N/2,N/2);
U(N/2+1:N,1:N/2,N/2+1:N) = ones(N/2,N/2,N/2);
U(N/2+1:N,N/2+1:N,1:N/2) = ones(N/2,N/2,N/2);

% %%%%%%%%%%%%%%%%%%% hyperbolique nappe

elseif (no_cond==15)
X1 = (X1+0*0.00).^2;
X2 = X2.^2;
X3 = X3.^2;    
R = (0.39)/(1+1/sqrt(2));
%R = 0.25;
R1 = R;
%R0 = 0.15;
R0 = R/sqrt(2);
U = sqrt((sqrt(X1 +X2 )-R1).^2 + X3) - R0;
U = 1/2 - 1/2*tanh(U/epsilon);



elseif (no_cond==16)
X1 = (X1+0*0.00).^2;
X2 = X2.^2;
X3 = X3.^2;
R0 = 0.1;
R1 = 0.39-R0;
U = sqrt((sqrt(X1 +X2 )-R1).^2 + X3) - R0;
U = 1/2 - 1/2*tanh(U/epsilon/2);


elseif (no_cond==17)
X1 = (X1+0*0.00);
X2 = X2;
R = (0.39)/(1+1/sqrt(2));
R1 = R
%R0 = 0.15;
R0 = R/sqrt(2);
U = sqrt((X1).^2 +X2.^2 )-R0/2;
U = 1/2 - 1/2*tanh(U/epsilon/2);

elseif (no_cond==18)
U(N/8+1:7*N/8,N/4+1:3*N/4,N/2-N/16+1:N/2+N/16 )= ones(3*N/4,N/2,N/8);
U(2*N/8+1:3*N/8,N/2-N/16+1:N/2+N/16,N/2-N/16+1:N/2+N/16 ) = zeros(N/8,N/8,N/8);
U(5*N/8+1:6*N/8,N/2-N/16+1:N/2+N/16,N/2-N/16+1:N/2+N/16 ) = zeros(N/8,N/8,N/8);

elseif (no_cond==19)
U(N/2-N/8+1:N/2+N/8,N/2-N/8+1:N/2+N/8,N/4+1:3*N/4)= ones(N/4,N/4,N/2);
U(1:N,N/2-N/16+1:N/2+N/16,N/2-N/16+1:N/2+N/16) = zeros(N,N/8,N/8);
U(N/2-N/16+1:N/2+N/16,1:N,N/2-N/16+1:N/2+N/16) = zeros(N/8,N,N/8);


elseif (no_cond==20), %%%%%%%%%%%% sphere trouée %%%%%%%%%
R0 = 0.43;
R1 = 0.43;
temp= max(sqrt(X1.^2+X2.^2+X3.^2)-R0, R1 - sqrt((X1-0.05).^2+X2.^2+(X3-0.2).^2));
U = 1/2 - 1/2*tanh(temp/(2*epsilon));
    

elseif (no_cond==21), 
R0 = 0.3;

temp= sqrt(X1.^2+X2.^2+3*X3.^2)-R0;
U = 1/2 - 1/2*tanh(temp/(2*epsilon));


elseif(no_cond ==22),
    
  l1 = 0.3;
  l2 = 0.2;
  l3 = 0.1;
  temp = max(max(max(max(max(X1 - l1,-X1-l1),X2-l2),-X2-l2),X3-l3),-X3-l3);
  U = 1/2 - 1/2*tanh(temp/(2*epsilon));
   

elseif(no_cond ==23),
   
  l3 = 0.35;
  R = 0.1;
  temp = max(max(sqrt(X1.^2 + X3.^2) - R,X2-l3),-X2-l3);
  U = 1/2 - 1/2*tanh(temp/(2*epsilon));
  %R = 0.4;
  %l3 = 0.125;
  %temp = max(max(sqrt(X1.^2 + 1.5*X2.^2) - R,X3-l3),-X3-l3);
  
  %%%%%%%%%%%%%%%%%%%%%% dumblelle
elseif (no_cond==25)
R_1 = 0.10;
R_2 = 0.3;
R_3 = 0.15;

U_temp = max((sqrt(X1.^2 +X2.^2 )-R_1), X3.^2 - R_2^2);
U_temp = min(min(U_temp,sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3-R_2).^2)-R_3),sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3+R_2).^2)-R_3);

U = 1/2 - 1/2*tanh(U_temp/epsilon/2);%U = 1/2 - 1/2*tanh(temp/(2*epsilon));
  
  
  
  
else
X1 = (X1).^2;
X2 = X2.^2;
X3 = X3.^2;  


R1=1;
R0 = 0.2;
d = max(X2/R1 + X3/R1 - X1-R0^2,  X1 - 0.3^2);
U1 = 1/2 - 1/2*tanh(d/(epsilon/2));
%U = zeros(N,N,N);
%U(N/4+1:N-N/4 ,N/4+1:N-N/4,N/4+1:N-N/4) = U1(N/4+1:N-N/4 ,N/4+1:N-N/4,N/4+1:N-N/4);
U = U1;
end




