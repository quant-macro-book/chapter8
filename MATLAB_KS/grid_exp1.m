function [grid] = grid_exp1( xmin, xmax, num_grid )
% Function grid_exp1
%  [grid] = grid_exp1( xmink, xmax, num_grid )
%
% Purpose:
%  Generate exponentially-spaced grids.
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  02/20/2013  T. Yamada   Translation: original in Carroll (2012)

    zmax = log(xmax+1.0);
    mesh = linspace( xmin, zmax, num_grid);
    grid = exp(mesh)-1.0;

return;