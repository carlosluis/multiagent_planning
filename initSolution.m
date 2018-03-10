function [p] = initSolution(po,pf,h,K)

diff = pf - po;
t = 0:h:(K-1)*h;

for i = 1:length(t)
    p(:,i) = po + t(i)*diff;
end
end
