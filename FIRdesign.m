% h = FIRdesign(wl, wu, N)

function [h, n] = FIRdesign(wl, wu, N)

R = rem(N,2);

if R ~= 0
    n = -(N-1)/2:(N-1)/2;
else
    n = -N/2:N/2;
end

h = (wu/pi)*sinc(wu*n/pi)-(wl/pi)*sinc(wl*n/pi);

end