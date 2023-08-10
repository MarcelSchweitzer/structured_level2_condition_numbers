%SURFACE_TEST_CORTH     Trying to get surface.

randn('state', 7)
n = 18; % number of matrix sizes tested
y = zeros(1,n); % for storing matrix sizes
y(1:10) = linspace(10,100,10);
y(11:15) = linspace(120,200,5);
y(16:18) = linspace(300,500,3);

num1 = 40; % max no. of G-reflectors applied
result_numbs = zeros(n,num1+1); result_poly = zeros(n,num1+1);

for k = 1:n
    
    y(k)
    [result_numbs(k,:),x] = test_corth(y(k),num1);
    % ^ row k stores condition numbers for matrix size y(k)
    
    p = polyfit(x,log(result_numbs(k,:)),4);
    result_poly(k,:) = polyval(p,x);

end

[X,Y] = meshgrid(x,y);
subplot(121)
mesh(x,y,log(result_numbs));
xlabel('No. of G-reflectors'), ylabel('Size of matrix (n)')
zlabel('log(cond(A))')
subplot(122)
mesh(x,y,result_poly);
xlabel('No. of G-reflectors'), ylabel('Size of matrix (n)')
zlabel('log(cond(A))')
title('Average condition number vs matrix size and no. of G-reflectors.')

fid = fopen('output_result_numbs','w');
fprintf(fid, '%g\n', result_numbs); % saves data by column
fclose(fid);
fid = fopen('output_result_poly','w');
fprintf(fid, '%g\n', result_poly);
fclose(fid);
fid = fopen('output_x','w');
fprintf(fid, '%g\n', x);
fclose(fid);
fid = fopen('output_y','w');
fprintf(fid, '%g\n', y);
fclose(fid);

% To get data back, use (for matrices):
% fid = fopen('output_result_numbs','r');
% A = fscanf(fid,'%g');
% fclose(fid);
% A = reshape(A,n,num1+1);
% OR (for vectors)
% fid = fopen('output_x','r');
% a = fscanf(fid,'%g');
% fclose(fid);
% a = a';
