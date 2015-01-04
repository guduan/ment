%--------------------------------------------------------
%
%    Matlab function to reconstruct image from sinogram
%
%       using the Convolution Back Projection method.
%
%    Kai Hock
%    Cockcroft Institute
%    10 June 2010
%
%--------------------------------------------------------------------
%
%   function recon_FBP(x0, y0, sinogram, angles, positions, centre, weights, nproj, npos)
%
%   Input:  x0          - x positions for reconstruction
%           y0          - y positions for reconstruction
%
%   Output: ReconFile   - name of file for reconstructed image
%
%
%--------------------------------------------------------------------
%
%
% nproj     - number of projections;
% angles    - column vector [nproj x 1] of angles
% weights   - column vector [nproj x 1] of weights for each projection.
%             For linear interpolation, use the sum of the angular intervals
%              on the two sides of each projection as the weight.
%             (Set to angle intervals.)
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

function recon_FBP(x0, y0, sinogram, angles, positions, centre, weights, nproj, npos)

% prange = max(positions', 'r') - min(positions', 'r');   %  vector [1 x nproj], range of positions
prange = max(positions', [],1) - min(positions', [],1);   %  vector [1 x nproj], range of positions

n1 = ceil(log(npos)/log(2));
N0 = 2^n1;                % smallest power of 2 greater than npos
N = 2^(n1+1);             % more zero padding

dr = prange / (npos-1);    % vector [1 x nproj] interval in position

rmax = dr * N;             % vector [1 x nproj] range of r

%-----FFT-------------------------------------------------

fprintf('FFT ... \n');
dw   = 1 ./ rmax;          % vector [1 x nproj]

sinogram1 = zeros(N, nproj);         %  zero padding
sinogram1(1:npos, :) = sinogram';   % each column is projections for one angle

%fpdata = fft(sinogram1); %  FFT every column

fpdata = zeros(N, nproj);

for iproj = 1:nproj;
    fpdata(:, iproj) = fft(sinogram1(:, iproj));
end



%----filter---------------------------------------
%  impulse response of bandlimited filter
%      Kak and Slaney (1988) p. 72
h0 = zeros(N, 1);

% compute impulse response without 1/tau^2 factor
%
% (note: tau is position interval dr,
%        which can be different for different angles)
%
% (N0 = N/2 is used instead of N for zero padding of h0
%  to avoid interperiod interference)
%
h0(1) = 1/4;                                     % n = 0
%h0(2:2:(N/2)) = -1 ./(((2:2:(N/2)) - 1) .^2 * pi^2);   %  n odd
for i = 2:2:(N/2),
    h0(i) = -1 /((i - 1)^2 * pi^2);   %  n odd
end

h0(N:-2:(N-N/2+2)) = h0(2:2:(N/2));


% create nproj columns of h0 and
% multiply each row by a corresponding 1/tau^2
%
% then multiply by tau again to get the right overall factor
% for the combined Fourier transforms
%
H1 = fft(h0);

% FFT ...
H2 = H1 * (ones(1, nproj) ./ dr);

% and filter
ffpdata = H2 .* fpdata;

%----IFFT-----------------------------------------

fprintf('IFFT ... \n');
%fdata1 = real(ifft(ffpdata))';   % ifft every column

fdata1 = zeros(nproj, N);

for iproj = 1:nproj;
    fdata1(iproj, :) = real(ifft(ffpdata(:, iproj)))';
end

%--back projection --------------------------------

fprintf('Back projection ... \n');

[yy, xx] = ndgrid(-y0, x0);    % for reconstruction
[nx, ny] = size(xx);
fdata = zeros(nx, ny);

for iproj = 1:nproj,
    fprintf('   angle %d \n', iproj);
    theta = angles(iproj);
    
    rho1 = xx*cos(theta) + yy*sin(theta);
    r2 = positions(iproj, :) - centre(iproj);       %  correct centre error
    
    fdata2 = fdata1(iproj, 1:npos);
    fdata4 = interp1(r2, fdata2, rho1, 'linear', 0);
    
    %  [m, n] = size(fdata)
    %  [m1, n1] = size(fdata4)
    %  mprintf('   %d  %d \n', m, n);
    %  mprintf('   %d  %d \n', m1, n1);
    
    fdata = fdata + fdata4 * weights(iproj);   % linear interpolation over angle
end

%--save reconstruction data ---------------

f0 = fdata/max(max(fdata));

u=fopen('recon_FBP.txt','w');
for i = 1:nx,
    for j = 1:ny,
        fprintf(u,'%6.3f   %6.3f   %6.3f\n', xx(i, j), yy(i, j), f0(i, j));
    end
    fprintf(u,'\n');
end
fclose(u);


fprintf('FBP OK!\n')
end
