function [p,v] = propState(po, a,A_p,A_v,K)

vo = [0 0 0];

new_p = A_p*a;
new_v = A_v*a;

p = [po'; new_p + repmat(po',K-1,1)];
v = [vo'; new_v];
end