
sinogram2=sinogram;

% FBP
filetext = 'projections0.txt';
isample = 1;
sino_txt(filetext, nproj, angles, weights, npos, positions, centre, sinogram2, isample);
% sino_txt(filetext, nproj, angles, weights, npos, positions, centre, sinogram, isample);
x0 = rho1:.01:rho2;
y0 = x0;
% recon_FBP(x0, y0, sinogram, angles, positions, centre, weights, nproj, npos);
recon_FBP(x0, y0, sinogram2, angles, positions, centre, weights, nproj, npos);

% show Projections.
projsdata=importdata('projections0.txt');
figure;
subplot(221);
plot(projsdata(:,1),projsdata(:,2));
title('Projections');
xlabel('x');ylabel('Projections');

%show FBP
fbpdata=importdata('recon_FBP.txt');
length_recon_FBP=sqrt(length(fbpdata));
fbp=reshape(fbpdata(:,3),length_recon_FBP,length_recon_FBP)';

subplot(222);
imshow(fbp,'XData',-1:0.01:1,'YData',1:-0.01:-1);axis on;colorbar;
title('ReCons\_FBP');set(gca,'ydir','normal')
xlabel('x');ylabel('y');

% MENT
ProjectFile = 'sinogram.bin';
fwrite_sinogram_1(ProjectFile, nproj, angles, weights, npos, positions, centre, sinogram2);
 

dos('./ment4c.apple 2 200');

mentdata=importdata('recon_MENT.txt');% file get from ment4c
length_recon_MENT=sqrt(length(mentdata));
ment=reshape(mentdata(:,3),length_recon_MENT,length_recon_MENT);

subplot(223);
imshow(ment,'XData',-1:0.01:1,'YData',1:-0.01:-1);axis on;colorbar;
title('ReCons\_MENT');set(gca,'ydir','normal')
xlabel('x');ylabel('y');
