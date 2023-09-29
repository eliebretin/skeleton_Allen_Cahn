N= 2^10;
x = linspace(-1,1,N);
[X1,X2] = meshgrid(x,x);

k = [0:N/2,-N/2+1:-1];
[K1,K2] = meshgrid(k,k);


I_aff = [N/16:N-N/16];/im

sigma_vec = [0.005,0.01,0.02];
for E_exemple =1:3,


%%%%%%%%%%%% 1exemple ,  tore %%%%%%%%%%%%%
if (E_exemple == 1)

r0 = 0.3;
r1 = 0.7;
dist_E = max(sqrt(X1.^2 + X2.^2) - r1, r0 - sqrt(X1.^2 + X2.^2));


elseif (E_exemple == 2)

dist_E = max(X1 -0.3,  -0.3 - X1);

elseif (E_exemple == 3)

dist_E = max(max(max(X1 -0.4,  -0.4 - X1),X2 - 0.15),-0.15 - X2);



end

E = dist_E <= 0;

clf;
imagesc(x(I_aff),x(I_aff),E(I_aff,I_aff));
title(['E']);
colorbar;
name_fig = ['E_',num2str(E_exemple),'.eps'];
print('-depsc', name_fig)


clf;
imagesc(x(I_aff),x(I_aff),dist_E(I_aff,I_aff));
title(['E']);
colorbar;
name_fig = ['distE_',num2str(E_exemple),'.eps'];
print('-depsc', name_fig)


for i_sigma = 1:3,


%%%%%%%%%%% calcul de S %%%%%%%%%%%%%%%%%%%

[n1,n2] = gradient(dist_E,2/N,2/N);
sigma = sigma_vec(i_sigma);
kernel = exp(-pi*sigma^2*(K1.^2 + K2.^2));
n1_sigma = ifft2(kernel.*fft2(n1));
n2_sigma = ifft2(kernel.*fft2(n2));

[partial1_n1,partial2_n1] = gradient(n1_sigma,2/N,2/N);
[partial1_n2,partial2_n2] = gradient(n2_sigma,2/N,2/N);

S_nsigma = partial1_n1.*n1_sigma.*n1_sigma + partial1_n2.*n1_sigma.*n2_sigma + partial2_n1.*n1_sigma.*n2_sigma + partial2_n2.*n2_sigma.*n2_sigma;


J_sigma = partial1_n1.* partial2_n2 - partial1_n2.*partial2_n1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
imagesc(x(I_aff),x(I_aff),S_nsigma(I_aff,I_aff));
title(['Sn_{\sigma} with \sigma = ',num2str(sigma)]);
colorbar;
name_fig = ['Sn_E_',num2str(E_exemple),'_sigma_',num2str(i_sigma),'.eps'];
print('-depsc', name_fig)




clf;
imagesc(x(I_aff),x(I_aff),J_sigma(I_aff,I_aff));
title(['Jn_{\sigma} with \sigma = ',num2str(sigma)]);
colorbar;
name_fig = ['Jn_E_',num2str(E_exemple),'_sigma_',num2str(i_sigma),'.eps'];
print('-depsc', name_fig)





end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
