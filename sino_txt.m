function sino_txt(filetext, nproj, angles, weights, npos, positions, centre, sinogram, isample)

fprintf('sino_txt\n');

u=fopen(filetext,'w');
for i = 1:isample:nproj
  for j = 1:npos
    rho = positions(i, j) - centre(i);
    fprintf(u,'%6.3f   %6.3f\n', rho, sinogram(i, j));
  end
  fprintf(u,'\n');
end
fclose(u);

 
end
