function [r,k] = reverse(x,n)

% x is the original sequence
% n is the time index for original sequence
% r is the reversed sequence
% k is the time index for the reversed sequence

p = size(x,2); %p is the total size(column) of the sequence
r = x(p:-1:1); %reverses the sequence from the last value
k= -n(p:-1:1); %reverses the time index

end