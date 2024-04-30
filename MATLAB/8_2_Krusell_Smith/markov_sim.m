function seq_state = markov_sim(tran, ng, ns)
%  [k_path, z_path] = law_of_motion( )
%
% Purpose:
%  Simulated path of Markov process.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 11/01/2019  T. Yamada   Original code

%% set random_seed
rng(225);
X = rand(1, ns);

%% set initial state
seq_state = zeros(ns, 1);
seq_state(1) = 1;
current = seq_state(1);

% main simulation part
for j = 2:ns
    next_pr = tran(current, 1);
    for i = 1:ng
        if X(j) <= next_pr
            seq_state(j) = i;
            if seq_state(j) ~= 0
                break
            end
        else
            next_pr = next_pr + tran(current,i+1);
        end
    end
    current = seq_state(j);
end

return;
