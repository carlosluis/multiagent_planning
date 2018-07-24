function [p,v] = propStatedmpc(po, vo,  a, A_initp, A_p, A_v)
K = length(a)/3;
new_p = A_p*(a) + A_initp*([po';vo']);
new_v = A_v*(a) + repmat(vo',K,1);

p = new_p;
v = new_v;
end