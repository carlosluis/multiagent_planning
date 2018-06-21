function inbounds = is_inbounds(p,pmin,pmax)
tol = 50e-3;
up = max(p(1,:)) < pmax(1)+tol && max(p(2,:)) < pmax(2)+tol && max(p(3,:)) < pmax(3)+tol;
down = min(p(1,:)) > pmin(1)-tol && min(p(2,:)) > pmin(2)-tol && min(p(3,:)) > pmin(3)-tol; 
inbounds = up && down;
end