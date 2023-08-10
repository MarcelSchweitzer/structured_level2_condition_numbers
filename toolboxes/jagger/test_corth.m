function [c,count] = test_corth(n,num1)
%TEST_CORTH    Testing function rand_corth.m BETTER
%   For a particular matrix size n:
%   count is the vector of number of G-reflectors applied
%   c is the vector of average condition numbers
%   num1 is the max number of G-reflectors applied
%   num2 is the number of tests at each level (then averaged)

num2 = 5;
c = zeros(1,num1+1); count = c;

for i = 1:num1+1
    for k = 1:num2
        A = rand_corth(n,i-1);
        c(i) = c(i) + cond(A);
    end
    count(i) = i-1;
end
c = c/num2;