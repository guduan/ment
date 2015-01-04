function fwrite_sinogram_1(filename, nproj, angles, weights, npos, positions, centre, sinogram) 
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

% fid = mtlb_fopen(filename, 'wb');
% 
% mtlb_fwrite(fid, nproj, 'int');
% mtlb_fwrite(fid, npos, 'int');
% 
% mtlb_fwrite(fid, angles, 'double');
% mtlb_fwrite(fid, weights, 'double');
% mtlb_fwrite(fid, centre, 'double');
% 
% mtlb_fwrite(fid, positions', 'double');
% mtlb_fwrite(fid, sinogram', 'double');
% 
% mclose(fid);

fid = fopen(filename, 'w');

fwrite(fid, nproj, 'int');
fwrite(fid, npos, 'int');

fwrite(fid, angles, 'double');
fwrite(fid, weights, 'double');
fwrite(fid, centre, 'double');

fwrite(fid, positions', 'double');
fwrite(fid, sinogram', 'double');

fclose(fid);

disp('sinogram has saved to file!');
end