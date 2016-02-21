function [w, n] = RectangularWindow2(N)

R = rem(N,2);

if R ~= 0
    n = -(N-1)/2:(N-1)/2;
else
    n = -N/2:N/2;
end

w = ones(1,N);

end