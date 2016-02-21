%{
    Unitstep function: u = unitstep(n)
    The unitstep function behaves like the unitstep function.

    input(n): time index

    output(u): unitstep sequence
%}

function u = unitstep(n)

for i = 1:1:length(n)
    if n(i) < 0
        u(i) = 0;
    else
        u(i) = 1;
    end
    
end