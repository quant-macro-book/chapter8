function [grid] = grid_exp3( xmin, xmax, num_grid )
% Function grid_exp3
%  [grid] = grid_exp3( xmin, xmax, num_grid )
%
% Purpose:
%  Generate triple exponentially-spaced grids.
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  02/20/2013  T. Yamada   Translation: original in Carroll (2012)

    zmax = log(log(log(xmax+1)+1)+1);
    mesh = linspace( xmin, zmax, num_grid);
    grid = exp(exp(exp(mesh)-1.0)-1.0)-1.0;

return;