function violation = CheckforColl (p,l,k,E1,r_min,order)
violation = false;
if (~isempty(l))
    for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist = norm(E1*(p-pj(:,k)),order);
        if dist < r_min
            violation = true;
            break;
        end
    end  
end
end