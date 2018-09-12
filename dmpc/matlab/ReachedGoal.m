function pass = ReachedGoal(p,pf,length_t,error_tol,N)

if(N>1)
    differ = squeeze(p(:,length_t,:))- squeeze(pf);
    max_dist = max(sqrt(sum(differ.^2,1)));
else
    differ = p(:,length_t) - pf';
    max_dist = max(sqrt(sum(differ.^2,1)));
end
    
pass = max_dist < error_tol;