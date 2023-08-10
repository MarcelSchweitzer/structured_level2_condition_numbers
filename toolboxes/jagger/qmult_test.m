%QMULT_TEST     Comparing speed of QMULT methods
%   Scores for Method 0 (Householder) in s0
%   Scores for Method 1 (QR) in s0

x = zeros(1,19); s0 = x; s1 = x;
x(1:9) = linspace(10,90,9);
x(10:19) = linspace(100,1000,10);
randn('state', 1)

for i = 1:19
    
    t = cputime;
    A = gallery('qmult',x(i),0);
%    A = qmult_unit(x(i),0);
    s0(i) = cputime - t;
    s = cputime;
    A = gallery('qmult',x(i),1);
%    A = qmult_unit(x(i),1);
    s1(i) = cputime - s;
    i
    
end  

warning off MATLAB:divideByZero
s2 = s0./s1;
plot(x,s2)
title('Comparing speeds of QMULT methods')
%title('Comparing speeds of QMULT\_UNIT methods')
xlabel('Matrix size (n)')
ylabel('Householder method time/QR method time')

%For QMULT:
%s2 = [ 4.1337 5.7726 5.9932 5.9110 5.5484 5.6241 5.4375 ...
%       5.2966 5.6596 6.2063 10.9953 15.5303 18.2030 20.4194  22.4324   ...
%       23.5556 26.4324 27.2765 27.7194];

%For QMULT_UNIT:
%s2 = [ 11.5652 12.2688 12.2979 9.9085 8.6594 8.3466 7.4652 7.7180 7.6686 ...
%       8.2263 9.0790 11.3179 14.2538 16.3985 17.9056 21.0870 46.6810 ...
%       33.9277 45.7881 ];