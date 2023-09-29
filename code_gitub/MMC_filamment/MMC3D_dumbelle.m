clear;

J=7;
N = [2^J,2^J,2^J];
L = [1,1,1];
x = linspace(-1/2,1/2,N(1));
 
 [X1,X2,X3] = meshgrid(x,x,x);
 
 epsilon = 2/N(1);
 
 bool_penalisation = 0;
 
 U1 = condition_initiale3d(N(1),epsilon,25);


 vv = VideoWriter(['test_dumbbell_',num2str(bool_penalisation),'.avi'],'Motion JPEG AVI');
    vv.Quality = 90;
    open(vv);



affiche_solution_3d2(x,U1,0*U1);
view(200,20)
pause(1);

 
 vol = sum(U1(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 dt=  epsilon^2/4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1_plus = [2:N(1),1];
i1_moins = [N(1),1:N(1)-1];
i2_plus = [2:N(2),1];
i2_moins = [N(2),1:N(2)-1];
i3_plus = [2:N(3),1];
i3_moins = [N(3),1:N(3)-1];

h1 = L(1)/N(1);
h2 = L(2)/N(2);
h3 = L(3)/N(3);



k = [0:N(1)/2,-N(1)/2+1:-1];
[K1,K2,K3] = meshgrid(k,k,k);
Lap = -(4*pi^2*(abs(K1).^2 + abs(K2).^2 +abs(K3).^2 ));
W_prim  = @(U) U.*(U-1).*(2*U-1);
W_seconde  = @(U) 6*U.*(U-1) + 1;
sqrtU = @(U) abs(U.*(1-U));
 alpha = 2;
 
 
 L = Lap - alpha/epsilon^2;
 Op_N  = @(U)  (-W_prim(U) + alpha*U)/epsilon^2;
 M = 1./(1 - dt*L);
sigma = 0.1/N(1);
Kernel_sigma = exp(-4*pi^2*sigma^2*(K1.^2 + K2.^2 + K3.^2));

if   bool_penalisation ==1; 
T = 0.1;
else
    T = 0.02
end

T_vec = linspace(0,0.99*T,10);
T_vec = [T_vec,2*T];
j_sauvegarde = 1;


for i=1:T/dt,
i
temps(i) = i*dt;
U = U1;

U = ifftn(Kernel_sigma.*fftn(U));


if bool_penalisation == 1,
%%%%%%%%%%%%%%%%%% div
gradU_1 = (U(i1_plus,:,:) - U(i1_moins,:,:))/(2*h1);
gradU_2 = (U(:,i2_plus,:) - U(:,i2_moins,:))/(2*h2);
gradU_3 = (U(:,:,i3_plus) - U(:,:,i3_moins))/(2*h3);

norm_grad = sqrt(gradU_1.^2 + gradU_2.^2 + gradU_3.^2 + 10^(-16));
A1 = gradU_1./norm_grad;
A2 = gradU_2./norm_grad;
A3 = gradU_3./norm_grad;   

%%%%%%%%%%%%%%%% regularisation avec sigma %%%%%%%%%%ùù

A1 = ifftn(Kernel_sigma.*fftn(A1));
A2 = ifftn(Kernel_sigma.*fftn(A2));
A3 = ifftn(Kernel_sigma.*fftn(A3));

gradU_21 = (A1(:,i2_plus,:) -  A1(:,i2_moins,:))/(2*h2);
gradU_11 = (A1(i1_plus,:,:) - A1(i1_moins,:,:))/(2*h1);
gradU_31 = (A1(:,:,i3_plus) -  A1(:,:,i3_moins))/(2*h3);

gradU_22 = (A2(:,i2_plus,:) -  A2(:,i2_moins,:))/(2*h2);
gradU_12= (A2(i1_plus,:,:) - A2(i1_moins,:,:))/(2*h1);
gradU_32 = (A2(:,:,i3_plus) -  A2(:,:,i3_moins))/(2*h3);


gradU_23 = (A3(:,i2_plus,:) -  A3(:,i2_moins,:))/(2*h2);
gradU_13= (A3(i1_plus,:,:) - A3(i1_moins,:,:))/(2*h1);
gradU_33 = (A3(:,:,i3_plus) -  A3(:,:,i3_moins))/(2*h3);


temp =  gradU_11.*A1.*A1 + gradU_21.*A2.*A1 +  gradU_31.*A3.*A1 ... 
     +  gradU_12.*A1.*A2 + gradU_22.*A2.*A2 +  gradU_32.*A3.*A2 ...
     +  gradU_13.*A1.*A3 + gradU_23.*A2.*A3 +  gradU_33.*A3.*A3 ;  

temp =  abs(16*(U.*(1-U)).^2).*ifftn((exp(+0.1*epsilon^2*Lap).*fftn(abs(temp))));


coef = 0.75/2*epsilon*N(1)^3;

temp1 = coef*temp + 1/epsilon^2;

else
 temp1 =  1/epsilon^2;
end


max_pen = max(temp1(:));

Op_N  = @(U)  - (W_prim(U).*temp1 - max_pen*U) ;
L = Lap  - max_pen;
M = 1./(1 - dt*L);

U1 =ifftn(M.*(fftn(U1 + dt*Op_N(U1))));

%volU1(i) = sum(U1(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (mod(i,50)==1)
    clf;
  affiche_solution_3d2(x,U1,0*U1);
  axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]) 
  view(200,20)
   
 frame = getframe(gcf);
 writeVideo(vv,frame);
pause(0.1);    
    
end


   
    if (i*dt > T_vec(j_sauvegarde))
       
           clf;
  affiche_solution_3d2(x,U1,0*U1);
  axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]) 
  view(200,20)
       name_title = ['t = ',num2str(i*dt)];
       title(name_title,'linewidth',2)
       axis square;
       
       name_fig = ['Test_dumbbell_',num2str(bool_penalisation),'_',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;


       
    end





end

close(vv)

