% clear all;clc;
rho1  = -1;     % lower limit of position in each projection
rho2  = 1;      % upper limit of position in each projection

nproj = 6;    % number of projections (or angles)
npos  = 100;    % number of rays (or positions)

[sinogram, angles, positions, centre, weights] = phantom_2b(rho1, rho2, nproj, npos) ;
% a={'ro','bo','go','ko','co','mo'};
a=['r','b','g','k','c','m'];

for n=1:6
    plot(sinogram(n,:),a(n));hold on;
end