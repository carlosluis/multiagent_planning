function [p,v,a] = initDMPC(po,pf,h,k_hor,K)
diff = pf-po;
K = K;
t = 0:h:(K-1)*h;

for i = 1:length(t)
    p(:,i) = po + 1*t(i)*diff/10;
end

p = p(:,1:k_hor);
v = zeros(3,k_hor);
a = zeros(3,k_hor);
end
