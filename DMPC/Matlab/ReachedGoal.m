function pass = ReachedGoal(p,pf,length_t,error_tol)

differ = squeeze(p(:,length_t,:))- squeeze(pf);
max_dist = max(sqrt(sum(differ.^2,1)));

pass = max_dist < error_tol;