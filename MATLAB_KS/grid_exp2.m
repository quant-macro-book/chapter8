function [grid] = grid_exp2( xmin, xmax, num_grid )
% Function grid_exp2
%  [grid] = grid_exp2( xmin, xmax, num_grid )
%
% Purpose:
%  Generate double exponentially-spaced grids.
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  02/20/2013  T. Yamada   Translation: original in Carroll (2012)

    zmax = log(log(xmax+1)+1);
    mesh = linspace( xmin, zmax, num_grid);
    grid = exp(exp(mesh)-1.0)-1.0;

return;