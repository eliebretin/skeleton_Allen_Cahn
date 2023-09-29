N= 2^7;
x = linspace(-1,1,N);
[X1,X2,X3] = meshgrid(x,x,x);

k = [0:N/2,-N/2+1:-1];
[K1,K2,K3] = meshgrid(k,k,k);


I_aff = [N/16:N-N/16];
sigma_vec = [0.005,0.01,0.02];


for E_exemple =1:3,
      

%%%%%%%%%%%% 1exemple ,  tore %%%%%%%%%%%%%
if (E_exemple == 1)

r0 = 0.3;
r1 = 0.7;
dist_E = max(sqrt(X1.^2 + X2.^2 + X3.^2) - r1, r0 - sqrt(X1.^2 + X2.^2 + X3.^2));

elseif (E_exemple == 2)
    
dist_E = max(X1 -0.3,  -0.3 - X1);

elseif (E_exemple == 3)
    
r0 = 0.1;
r1 = 0.7;
dist_E =  sqrt((sqrt(X1.^2 + X2.^2) - r1).^2 +  X3.^2) - 0.1;

    
end

E = dist_E <= 0;


%  clf
%  v = dist_E;
%  p = patch(isosurface(x,x,x,v,0));
%  isonormals(x,x,x,v,p)
%  set(p,'FaceColor','red','EdgeColor','none');
% 
%  daspect([1 1 1])
%  view(3); 
%  camlight
%  lighting gouraud
%  alpha (0.5)
% 
% title('E');
% name_fig = ['E_3D_',num2str(E_exemple),'.eps'];
% print('-depsc', name_fig)




for i_sigma = 3:3,


%%%%%%%%%%% calcul de S %%%%%%%%%%%%%%%%%%%

[n1,n2,n3] = gradient(dist_E,2/N,2/N,2/N);
%%% normalisation %%%%%ù
norm_n = sqrt(n1.^2 + n2.^2 + n3.^2);
n1 = n1./norm_n;
n2 = n2./norm_n;
n3 = n3./norm_n;



sigma = sigma_vec(i_sigma);

kernel = exp(-pi*sigma^2*(K1.^2 + K2.^2 + K3.^2));
n1_sigma = ifftn(kernel.*fftn(n1));
n2_sigma = ifftn(kernel.*fftn(n2));
n3_sigma = ifftn(kernel.*fftn(n3));



[partial1_n1,partial2_n1,partial3_n1] = gradient(n1_sigma,2/N,2/N,2/N);
[partial1_n2,partial2_n2,partial3_n2] = gradient(n2_sigma,2/N,2/N,2/N);
[partial1_n3,partial2_n3,partial3_n3] = gradient(n3_sigma,2/N,2/N,2/N);


S_nsigma = partial1_n1.*n1_sigma.*n1_sigma +partial1_n2.*n1_sigma.*n2_sigma + partial1_n3.*n1_sigma.*n3_sigma + ...
           partial1_n2.*n2_sigma.*n1_sigma +partial2_n2.*n2_sigma.*n2_sigma + partial2_n3.*n2_sigma.*n3_sigma + ...
           partial1_n3.*n3_sigma.*n1_sigma +partial3_n2.*n3_sigma.*n2_sigma + partial3_n3.*n3_sigma.*n3_sigma ; 
       
       

imagesc(x(I_aff),x(I_aff),S_nsigma(I_aff,I_aff,N/2));
title(['Sn_{\sigma} 3D on x_3 = 0 with \sigma = ',num2str(sigma)]);
colorbar;
name_fig = ['SnE_3D_x3_',num2str(E_exemple),'_sigma_',num2str(i_sigma),'.eps'];
print('-depsc', name_fig)


imagesc(x(I_aff),x(I_aff),squeeze(S_nsigma(I_aff,N/2,I_aff))');
title(['Sn_{\sigma} 3D on x_2 = 0 with \sigma = ',num2str(sigma)]);
colorbar;
name_fig = ['SnE_3D_x2_',num2str(E_exemple),'_sigma_',num2str(i_sigma),'.eps'];
print('-depsc', name_fig)

imagesc(x(I_aff),x(I_aff),squeeze(S_nsigma(N/2,I_aff,I_aff))');
title(['Sn_{\sigma} 3D on x_1 = 0 with \sigma = ',num2str(sigma)]);
colorbar;
name_fig = ['SnE_3D_x1_',num2str(E_exemple),'_sigma_',num2str(i_sigma),'.eps'];
print('-depsc', name_fig)



end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
