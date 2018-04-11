function violation = CheckCollEllipDMPC (p,l,n,k,E)
violation = false;
for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
    if(i~=n)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist = norm(E^(-1)*(p-pj(:,k)));
        violation = (violation) || (dist < 0.5);
    end
end  
end
