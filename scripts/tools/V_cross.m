function [c] = V_cross(a, b)
    % This is ~10 times quicker than the Matlab cross() function
    % for me (AMW). We assume the arguments are both 3-vectors and
    % avoid the checks and reshaping needed for the more general case.
    % (According to the prfiler, this moves cross from the most expensive
    % child function costing ~50% of the time to the third most expensive 
    % child costing ~10% of the time).

    c = [a(2).*b(3)-a(3).*b(2)
         a(3).*b(1)-a(1).*b(3)
         a(1).*b(2)-a(2).*b(1)];
    
return
