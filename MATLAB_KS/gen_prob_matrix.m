function [pi_dist, prob] = gen_prob_matrix( duration, unemp )
% Function gen_prob_matrix
%  [pi_dist, prob] = gen_prob_matrix( duration, unemp, ne )
%  
% Input:
%  duration: duration period of unemployment
%  unemp:    unemployment rate
%  ne:       #state (must be 2)
%
% Output:
%  pi_dist:  stationary distribution
%  prob:     transition probability matrix
%
%  Record of revisions:
%    Date      Programmer  Description of change
%  ==========  ==========  =====================
%  09/11/2019  T. Yamada   Original

% stationary distribution
pi_dist = zeros(2, 1);
pi_dist(1) = 1.0 - unemp;
pi_dist(2) = unemp;

% transition probability matrix
prob(2, 2) = (duration - 1.0) / duration;
prob(2, 1) = 1.0 - prob(2, 2);
prob(1, 1) = 1.0 - (1.0-prob(2, 2))*(pi_dist(2)/pi_dist(1));
prob(1, 2) = 1.0 - prob(1, 1);

return;