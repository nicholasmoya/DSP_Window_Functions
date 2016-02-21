function[y,n] = convolve(h,nh,x,nx)

    % h is the system's impulse response,
    % nh is the time index for h,
    % x is the input sequence,
    % nx is the time index for x,
    % y is the convolution of h and x and n is the time index for y.

    l = 1;

    for k = (min(nx)+min(nh)): (max(nx)+max(nh));

        [y1,n1] = reverse(h,nh);
        [y2,n2] = shift(y1,n1,k);
        [yn,nn] = compsigmul(x,nx,y2,n2);

        temp = 0;
        m = 1;
        a = 0;

        for i = min(nn):max(nn);
            
            a = temp + yn(m);
            temp = a;
            m = m + 1;
            
        end

        s(l) = a;
        l = l + 1;

    end

y = s;
n = (min(nx)+min(nh)): (max(nx)+max(nh));

end