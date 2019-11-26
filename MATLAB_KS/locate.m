function loc = locate(xx,x)
% Purpose:
% From Numerical Recipes.
%
% Record of revisions:
%    Date      Programmer   Description of change
% ==========   ==========   =====================
% 11/10/2010   T. Yamada    Original code

[r l] = size(xx);
ascnd = logical(xx(l) >= xx(1));
jl = 0;
ju = l + 1;
for i = 1:l
    if ju - jl <= 1
       break 
    end
    jm = floor((ju+jl)/2);
    if ascnd == logical(x>=xx(jm))
        jl = jm;
    else
        ju = jm;
    end
end
if x == xx(1)
    loc = 1;
elseif x == xx(l)
    loc = l - 1;
else
    loc = jl;
end
return;