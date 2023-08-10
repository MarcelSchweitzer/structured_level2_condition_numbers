%LOADMYOUTPUT   Load data.

n = 18; num1 = 40;

fid = fopen('output_result_numbs','r');
numbs = fscanf(fid,'%g');
fclose(fid);
numbs = reshape(numbs,n,num1+1);

fid = fopen('output_result_poly','r');
poly = fscanf(fid,'%g');
fclose(fid);
poly = reshape(poly,n,num1+1);

fid = fopen('output_x','r');
x = fscanf(fid,'%g');
fclose(fid);
x = x';

fid = fopen('output_y','r');
y = fscanf(fid,'%g');
fclose(fid);
y = y';