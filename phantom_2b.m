function [sinogram, angles, positions,centre,weights] = phantom_2b(rho1, rho2, nproj, npos) 

%--------------------------------------------
% deff('y = heaviside(x)',  ...  
%      '[m, n] = size(x); ...
%       y = zeros(m, n); ...
%       for i = 1:m, ...
%         for j = 1:n, ...
%           if x(i, j) > 0, ...
%             y(i, j) = 1; ...
%           end ...
%         end ...
%       end');
%---------------------------------------

%
%  To generate images and projections 
%  using ellipses for testing
%  tomography reconstruction codes
%
%  E.g. input:
%
%  rho1 = -1;        % minimum rho
%  rho2 = 1;         % maximum rho
%  nproj = 110;      % number of angles
%  npos  = 127;      % number of rays
%  PNG = 0;          % 1 save png, 0 don't save
%  EPS = 0;          % 1 save eps, 0 don't save
%  MAT = 0;          % 1 save mat, 0 don't save
%
%  phantom_2b('sino', -1, 1, 110, 127);
%
%  Kai Hock
%  Cockcroft Institute
%  10 June 2010
%---------------------------------------

%--describe  image -----
%
% Each row describes an ellipse.
% 'Refractive index' refers to intensity.
% At points where ellipses,overlap, the intensities are cumulated.
%
%-----------------------------------------------------------
%    center          major   minor   rotation   refractive
%    coordinate      axis    axis    angle      index
%     x      y                       (deg)

%  test 1: solid cylinder
AI1 = [0      0       0.7    0.7      0         2.0];

% test 2: hollow cylinder
AI2 = [0      0       0.7    0.7      0         2.0
       0      0       0.6    0.6      0        -2.0];
  
% test 2: hollow cylinder with rod
AI3 = [0      0       0.7    0.7      0         2.0
       0      0       0.6    0.6      0        -2.0
       0.2    0.2     0.2    0.2      0         2.0];

%AI = AI1;
%AI = AI2;
AI = AI3;

%---image size -----------

xlow  = -1; % image left
xhigh =  1; % image right
ylow  = -1; % image bottom
yhigh =  1; % image top
xstep = 0.01; % x resolution
ystep = 0.01; % y resolution

%------------------------------------------------------------

x1 = AI(:, 1);   % centre coordinate x
y1 = AI(:, 2);  % centre coordinate y
A  = AI(:, 3);   % semi major axis
B  = AI(:, 4);   % semi minor axis
a1 = AI(:, 5);   % (deg) rotation angle
ri = AI(:, 6);   % refractive index
n1 = length(x1); % number of ellipses

%----intensity matrix -------

x0 = xlow:xstep:xhigh;
y0 = ylow:ystep:yhigh;
nx = length(x0);
ny = length(y0);
[yy, xx] = ndgrid(-y0, x0);
f0 = zeros(ny, nx);

for i1 = 1:n1
%     mprintf("ellipse %d\n", i1);
  fprintf('ellipse %d\n', i1);
  alpha1 = a1(i1)/180*pi;
  x =  (xx-x1(i1))*cos(alpha1) + (yy-y1(i1))*sin(alpha1);
  y = -(xx-x1(i1))*sin(alpha1) + (yy-y1(i1))*cos(alpha1);
  f0 = f0 + ri(i1)*heaviside(1 - (x/A(i1)) .^2 - (y/B(i1)) .^2);
end
fmax = max(max(f0));
f0 = f0/fmax;

%--save phantom data ---------------

% u=mopen('phantom.txt','w');
u=fopen('phantom.txt','w');
for i = 1:nx,
  for j = 1:ny,
    fprintf(u,'%6.3f   %6.3f   %6.3f\n', xx(i, j), yy(i, j), f0(i, j));
  end
  fprintf(u,'\n');
end
fclose(u);


%--------------------------------------------------------------------
%
% ProjectFile data format:
%
% nproj     - number of projections;
% angles    - column vector [nproj x 1] of angles
% weights   - column vector [nproj x 1] of weights for each projection.
%             For linear interpolation, use the sum of the angular intervals
%              on the two sides of each projection as the weight.
%             (Set to ones by default.) 
% npos      - number of positions, must be the same for each projection
% positions - matrix [nproj x npos] for positions of each and every projection
%             (for each projection, the positions must have uniform intervals, 
%              and the same number of points)
% centre    - column vector [nproj x 1], 
%             the position of the centre of rotation for each projection
%             (set to zeros by default) 
% sinogram  - matrix of projection data [nproj x npos], 
%             where npos = number of positions across the projections
%
%--------------------------------------------------------------------

%-- set projection parameters ------------

centre  = zeros(nproj, 1);
weights = ones(nproj, 1)*pi/nproj;
angles = (1:nproj)'/nproj*pi;

drho = (rho2 - rho1)/(npos-1);
rho  = rho1 + (0:(npos-1))*drho;
for i = 1:nproj,
  positions(i, :) = rho;
end

%----- generate projections ------------------------
fprintf('Generating projections\n');

sinogram = zeros(nproj, npos);

for iproj = 1:nproj,
    p1 = zeros(1, npos);
    
    for i1 = 1:n1,
        theta = angles(iproj) - a1(i1)/180*pi;
        a2 = sqrt((A(i1)*cos(theta))^2 + (B(i1)*sin(theta))^2);
        
        s1 = sqrt(x1(i1)^2 + y1(i1)^2);
        gamma1 = atan2(y1(i1), x1(i1));
        
        for ipos = 1:npos,
            t1 = rho(ipos) - s1 * cos(gamma1 - angles(iproj));
            if (abs(t1) < a2),
                p1(ipos) = p1(ipos) + 2*ri(i1)*A(i1)*B(i1)/a2^2 * sqrt(a2^2 - t1^2);
            end
        end
    end
    sinogram(iproj, :) = p1;
end

%--save phantom data ---------------

sino = sinogram'/max(max(sinogram));

u=fopen('sinogram.txt','w');
for i = 1:npos,
  for j = 1:nproj,
    fprintf(u,'%6.3f   %6.3f   %6.3f\n', rho(i), angles(j)/pi*180, sino(i, j));
  end
  fprintf(u,'\n');
end
fclose(u);

fprintf('Phantom is OK!\n');

end
% endfunction