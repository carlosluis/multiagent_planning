function violation = CheckforAllColl (p,l,K,r_min,addconstr)
violation = false;
if (~isempty(l))
    for k = 1:K
        if (~ismember(k,addconstr))
            for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
                pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
                dist = norm(p(:,k)-pj(:,k));
                if dist < r_min 
                    violation = true;
                    return;
                end
            end
        end
    end
end
end