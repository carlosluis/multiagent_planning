function [p] = initAllSolutions(po,pf,h,K)
po = squeeze(po);
pf = squeeze(pf);
N = size(po,2);
t = 0:h:(K-1)*h;
for i = 1:N
   diff = pf(:,i) - po(:,i);
   for k = 1:length(t)
      p(:,k,i) =  po(:,i) + t(k)*diff/((K-1)*h);
   end
end