function violation = CheckforAllColl (p,l,K,E1,r_min,addconstr,order)
violation = false;
if (~isempty(l))
    for k = 1:K
        if (~ismember(k,addconstr))
            for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
                pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
                dist = norm(E1*(p(:,k)-pj(:,k)),order);
                if dist < r_min 
                    violation = true;
                    return;
                end
            end
        end
    end
end
end