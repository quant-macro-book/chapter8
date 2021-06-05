function [rent, wage] = factor_price(agg_cap, agg_lab, tfp, alpha, delta)
% Purpose:
% Assuming the production function is of Cobb=Douglas type,
% compute factor prices from the first order condition.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 11/09/2010  T. Yamada   Original code

rent = tfp*     alpha *agg_cap^(alpha-1.0)*agg_lab^(1.0-alpha) - delta;
wage = tfp*(1.0-alpha)*agg_cap^ alpha     *agg_lab^(   -alpha);

return;
