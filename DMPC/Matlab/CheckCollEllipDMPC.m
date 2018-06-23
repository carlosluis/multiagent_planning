function [violation,lim_eps] = CheckCollEllipDMPC (p,l,n,k,E1,rmin,order)
violation = false;
idx = 1;
for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
    if(i~=n)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist = norm(E1*(p-pj(:,k)),order);
        violation = (violation) || (dist < rmin);
        lim_eps(idx) = abs(rmin - dist)+ 0.05+k*(0.004)-0.004;
        idx = idx + 1;
    end  
end
end
