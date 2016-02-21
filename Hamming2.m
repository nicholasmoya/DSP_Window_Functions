

function [w, n] = Hamming2(N)

R = rem(N,2);

if R ~= 0
    n = -(N-1)/2:(N-1)/2;
    w = 0.54+0.46*cos(2*pi*n/(N-1));
else
    n = -N/2:N/2;
    w = 0.54+0.46*cos(2*pi*n/(N));
end

end