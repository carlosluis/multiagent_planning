function [p,v,a] = getStates(po,atot,A_p,A_v,K,N)

vo = [0;0;0];
po = squeeze(po);
for i = 1:N
   ai = atot((3*K*(i-1)+1):3*K*i);
   a(:,:,i) = vec2mat(ai',3)';
   new_p = A_p*ai;
   new_v = A_v*ai;
   pi = [po(:,i); new_p + repmat(po(:,i),K-1,1)];
   vi = [vo; new_v];
   p(:,:,i) = vec2mat(pi,3)';
   v(:,:,i) = vec2mat(vi,3)';
end