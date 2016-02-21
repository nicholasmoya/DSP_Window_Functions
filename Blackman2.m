

function [w, n] = Blackman2(N)

R = rem(N,2);

if R ~= 0
    n = -(N-1)/2:(N-1)/2;
    w = 0.42+0.5*cos(2*pi*n/(N-1))+0.08*cos(4*pi*n/(N-1));
else
    n = -N/2:N/2;
    w = 0.42+0.5*cos(2*pi*n/(N))+0.08*cos(4*pi*n/(N));
end

end