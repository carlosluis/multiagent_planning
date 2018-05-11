function violation = CheckCollDMPC (p,l,n,k,r_min)
violation = false;
for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
    if(i~=n)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist = norm(p-pj(:,k));
        violation = (violation) || (dist < r_min);
    end
end  
end
