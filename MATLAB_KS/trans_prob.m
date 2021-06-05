function [prob, tran_zz] = trans_prob( Params )
% Function trans_prob
%  [prob, tran_zz] = trans_prob( )
%  
% Input:
%  Params
%
% Output:
%  prob:    transition probability of (e,z)
%  tran_zz: transition probability of z
%
%  Record of revisions:
%    Date      Programmer  Description of change
%  ==========  ==========  =====================
%  09/11/2019  T. Yamada   Original

% transition probability of (z,z'): business cycle
[pi_z, tran_zz] = gen_prob_matrix(Params.durgd, 0.5);

% copy-n-paste from A.Smith's code
pgg00 = (Params.durug-1.0)/Params.durug; % tran_ee(2,2)|good
pbb00 = (Params.durub-1.0)/Params.durub; % tran_ee(2,2)|bad
pbg00 = 1.25*pbb00;
pgb00 = 0.75*pgg00;
pgg01 = (Params.unempg - Params.unempg*pgg00)/(1.0-Params.unempg);
pbb01 = (Params.unempb - Params.unempb*pbb00)/(1.0-Params.unempb);
pbg01 = (Params.unempb - Params.unempg*pbg00)/(1.0-Params.unempg);
pgb01 = (Params.unempg - Params.unempb*pgb00)/(1.0-Params.unempb);
pgg = (Params.durgd-1.0)/Params.durgd;       % tran_zz(1,1)
pgb = 1.0 - (Params.durbd-1.0)/Params.durbd; % tran_zz(1,2)

pgg10 = 1.0 - (Params.durug-1.0)/Params.durug; % tran_ee(1,2)|good
pbb10 = 1.0 - (Params.durub-1.0)/Params.durub; % tran_ee(1,2)|bad
pbg10 = 1.0 - 1.25*pbb00;
pgb10 = 1.0 - 0.75*pgg00;
pgg11 = 1.0 - (Params.unempg - Params.unempg*pgg00)/(1.0-Params.unempg);
pbb11 = 1.0 - (Params.unempb - Params.unempb*pbb00)/(1.0-Params.unempb);
pbg11 = 1.0 - (Params.unempb - Params.unempg*pbg00)/(1.0-Params.unempg);
pgb11 = 1.0 - (Params.unempg - Params.unempb*pgb00)/(1.0-Params.unempb);
pbg = 1.0 - (Params.durgd-1.0)/Params.durgd; % tran_zz(2,1)
pbb = (Params.durbd-1.0)/Params.durbd;       % tran_zz(2,2)

prob(1, 1, 1, 1) = pgg*pgg11;
prob(1, 1, 1, 2) = pbg*pbg11;
prob(1, 1, 2, 1) = pgg*pgg01;
prob(1, 1, 2, 2) = pbg*pbg01;
prob(1, 2, 1, 1) = pgb*pgb11;
prob(1, 2, 1, 2) = pbb*pbb11;
prob(1, 2, 2, 1) = pgb*pgb01;
prob(1, 2, 2, 2) = pbb*pbb01;
prob(2, 1, 1, 1) = pgg*pgg10;
prob(2, 1, 1, 2) = pbg*pbg10;
prob(2, 1, 2, 1) = pgg*pgg00;
prob(2, 1, 2, 2) = pbg*pbg00;
prob(2, 2, 1, 1) = pgb*pgb10;
prob(2, 2, 1, 2) = pbb*pbb10;
prob(2, 2, 2, 1) = pgb*pgb00;
prob(2, 2, 2, 2) = pbb*pbb00;

% debug: check consistency (=1)
% test11 = prob(1, 1, 1, 1)/tran_zz(1, 1) + prob(1, 1, 2, 1)/tran_zz(1, 1);
% test12 = prob(1, 1, 1, 2)/tran_zz(1, 2) + prob(1, 1, 2, 2)/tran_zz(1, 2);
% test21 = prob(1, 2, 1, 1)/tran_zz(2, 1) + prob(1, 2, 2, 1)/tran_zz(2, 1);
% test22 = prob(1, 2, 1, 2)/tran_zz(2, 2) + prob(1, 2, 2, 2)/tran_zz(2, 2);

return;
