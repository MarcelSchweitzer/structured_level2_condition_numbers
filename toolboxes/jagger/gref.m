%GREF   Investigating condition number of G-reflectors
%       For matrices of sizes 1:n - apply 1 G-reflector
%       results stores condition numbers

randn('state', 5)
n = 19;
y = zeros(n,1);
y(1:9) = linspace(10,90,9);
y(10:19) = linspace(100,1000,10);
c = zeros(1,n);

num1 = 30;  % how many to apply for average
num2 = 15;
num3 = 5;
results = zeros(9*num1+5*num2+5*num3,2);
count = 1;

for i = 1:9
    i
    for j = 1:num1
        A = rand_corth(y(i),1);
        results(count,1) = y(i);
        results(count,2) = cond(A);
        count = count + 1;
    end
end

for i = 10:14
    i
    for j = 1:num2
        A = rand_corth(y(i),1);
        results(count,1) = y(i);
        results(count,2) = cond(A);
        count = count + 1;
    end
end

for i = 15:19
    i
    for j = 1:num3
        A = rand_corth(y(i),1);
        results(count,1) = y(i);
        results(count,2) = cond(A);
        count = count + 1;
    end
end

A = zeros(9*num1+5*num2+5*num3,4);
A(:,1) = ones(9*num1+5*num2+5*num3,1);
A(:,2) = results(:,1);
A(:,3) = A(:,2).^2;
% [ 1 n n^2 ]
A(:,4) = log(results(:,2)); % i.e. b

X = lsqr(A(:,1:3),A(:,4),[],50); 
% gives LS solution [a0,a1,a2] s.t. 
% log(k2(A)) "=" a0 + a1*n + a2*n^2

plot(results(:,1),log(results(:,2)),'o')
xlabel('Size of matrix (n)'), ylabel('log(cond(G))')
title('Condition numbers of G-reflectors.')
hold
c = zeros(1,n); for k = 1:n, c(k) = X(1) + X(2)*y(k) + X(3)*y(k)^2; end
plot(y,c), hold off