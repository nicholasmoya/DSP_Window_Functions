function [s,n] = shift(x, N, nx)

%shift function
%x is the original sequence
%s is the shifted sequence
%N is the amount of shift
%nx is the time index for x

n = nx + N;
s = x;

end