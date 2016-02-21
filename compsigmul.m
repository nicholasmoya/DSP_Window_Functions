function [s, n] = compsigmul(x1,n1,x2,n2)

    % s1 and s2 are the equal length sequences corresponding to x1 and x2
    % n is the common time index

    n = min(min(n1),min(n2)): max(max(n1),max(n2));
    s1 = zeros(1,length(n));
    s2 = s1;

    s1(find((n >= min(n1)) & (n <= max(n1)) == 1)) = x1 ;% x1 with duration of s1
    s2(find((n >= min(n2)) & (n <= max(n2)) == 1)) = x2 ; % x2 with duration of s2

    s = s1.*s2;

end