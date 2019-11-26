function [next_cap] = approx_agg(agg_cap, state, icept, slope)
% Purpose:
% From approximate aggregation,
% compute the next period's aggregate capital.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 11/09/2010  T. Yamada   Original code

next_cap = exp(icept(state) + slope(state)*log(agg_cap));

return;