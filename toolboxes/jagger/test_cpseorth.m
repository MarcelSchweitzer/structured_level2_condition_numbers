%TEST_CPSEORTH    Testing RAND_CPSEORTH

randn('state', 7)
n = 20; % number of matrix sizes tested
y = zeros(n,1); % for storing matrix sizes
y(1:10) = linspace(10,100,10); 
y(11:15) = linspace(120,200,5);
y(16:18) = linspace(300,500,3);
y(19:20) = [750 1000];

num1 = 40; % max no. of G-reflectors applied
num3 = [1 5 10 15 20 25 30 35 40];
num2 = 5; % no. of trials at each n, k for n = 1:18
num4 = 3; % no. of trials at each n, k for n = 19:20

x = linspace(1,num1,num1); x = x';
results = zeros(18*num1*num2 + 18*num4,3);
% results = [ k ; n ; cond(A = rand_cpseorth(p,q,k)) ]
% with each num2 block of rows being same p,q,k,
% with p, q picked as q = floor(n/2); p = p-q;

count = 1;

for i = 1:18
    i
    for j = 1:num1
        for k = 1:num2
            
            q = floor(y(i)/2); p = y(i)-q;
            A = rand_cpseorth(p,q,j);
            results(count,1) = j;
            results(count,2) = y(i);
            results(count,3) = cond(A);
            count = count + 1;
        
        end
    end
end

for i = 19:20
    i
    for j = 1:9
        for k = 1:num4
            
            q = floor(y(i)/2); p = y(i)-q;
            A = rand_cpseorth(p,q,num3(j));
            results(count,1) = num3(j);
            results(count,2) = y(i);
            results(count,3) = cond(A);
            count = count + 1;
        
        end
    end    
end

results = sortrows(results,3);
count = 1;
while results(count,3) < 10E16
    count = count + 1;
end
count = count - 1;
count
results = results(1:count,:); % knock out #'s > 10E16

results = sortrows(results,1);
results = sortrows(results,2);

A = zeros(count,5);
A(:,1) = ones(count,1);
A(:,2:3) = results(:,1:2);
A(:,4) = A(:,2).^2;
% [ 1 k n k^2 ]
A(:,5) = log(results(:,3)); % i.e. b

fid = fopen('output_cpseorth_A','w');
fprintf(fid, '%g\n', A); % saves data by column
fclose(fid);

X = lsqr(A(:,1:4),A(:,5),[],50); 
% gives LS solution [a0,a1,a2,a3] s.t. 
% log(k2(A)) "=" a0 + a1*k + a2*n + a3*k^2

fid = fopen('output_cpseorth_LSsol','w');
fprintf(fid, '%g\n', X);
fclose(fid);

% now plots of n = 10;
%c = zeros(1,40); for k = 1:40, c(k) = X(1) + X(2)*k + X(3)*10 + X(4)*k^2; end
%plot(1:1:40,c), hold off

% To retreive data, use (for matrices):
% fid = fopen('output_cpseorth_A','r');
% B = fscanf(fid,'%g');
% fclose(fid);
% [m,n] = size(B);
% B = reshape(B,m/5,5);
% OR (for vectors)
% fid = fopen('output_cpseorth_LSsol','r');
% a = fscanf(fid,'%g');
% fclose(fid);