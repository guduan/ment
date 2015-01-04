%------------------------
%
%  Tomography test codes:
%
%  To generate phantom image, sinogram, and
%  to do reconstruction using the
%  Convolution Back Projection method.
%
%  Kai Hock
%  Cockcroft Institute
%  10 June 2010
%-----------------------
% Converted to Matlab codes by Duan Gu. 2014.03.18.
% Based on Kai's Scilab Codes.
%------------------------
clear all;clc;close all;


rho1  = -1;     % lower limit of position in each projection
rho2  = 1;      % upper limit of position in each projection

nproj = 6;    % number of projections (or angles)
npos  = 100;    % number of rays (or positions or num. of bins)

[sinogram, angles, positions, centre, weights] = phantom_2b(rho1, rho2, nproj, npos) ;

phantomdata=importdata('phantom.txt');
length_original=sqrt(length(phantomdata));
original=reshape(phantomdata(:,3),length_original,length_original)';

% show original x-y
figure;
subplot(221);
imshow(original,'XData',-1:0.01:1,'YData',1:-0.01:-1);axis on;colorbar;
title('Original');set(gca,'ydir','normal');
xlabel('x');ylabel('y');


% FBP
filetext = 'projections0.txt';
isample = 1;
sino_txt(filetext, nproj, angles, weights, npos, positions, centre, sinogram, isample);
x0 = rho1:.01:rho2;
y0 = x0;
recon_FBP(x0, y0, sinogram, angles, positions, centre, weights, nproj, npos);

% show Projections.
projsdata=importdata('projections0.txt');
subplot(222);
plot(projsdata(:,1),projsdata(:,2));
title('Projections');
xlabel('x');ylabel('Projections');

%show FBP
fbpdata=importdata('recon_FBP.txt');
length_recon_FBP=sqrt(length(fbpdata));
fbp=reshape(fbpdata(:,3),length_recon_FBP,length_recon_FBP)';


subplot(223);
imshow(fbp,'XData',-1:0.01:1,'YData',1:-0.01:-1);axis on;colorbar;
title('ReCons\_FBP');set(gca,'ydir','normal');
xlabel('x');ylabel('y');

% MENT
ProjectFile = 'sinogram.bin';
fwrite_sinogram_1(ProjectFile, nproj, angles, weights, npos, positions, centre, sinogram);


ostype=computer;
switch ostype
    case {'PCWIN','PCWIN64'}
        cmds='D:/ment/ment4c.exe 2 200';  % path to your ment4c.exe file. It should be re-compiled beforehand.
    case {'GLNX86','GLNXA64'}
        cmds='./ment4c.linux 2 200';
    case 'MACI64'
        cmds='./ment4c.apple 2 200';
end
dos(cmds);

mentdata=importdata('recon_MENT.txt');% file get from ment4c
length_recon_MENT=sqrt(length(mentdata));
ment=reshape(mentdata(:,3),length_recon_MENT,length_recon_MENT);

subplot(224);
imshow(ment,'XData',-1:0.01:1,'YData',1:-0.01:-1);axis on;colorbar;
title('ReCons\_MENT');set(gca,'ydir','normal');
xlabel('x');ylabel('y');
