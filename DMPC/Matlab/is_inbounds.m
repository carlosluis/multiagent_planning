function inbounds = is_inbounds(p,pmin,pmax)

up = max(p(1,:)) < pmax(1) && max(p(2,:)) < pmax(2) && max(p(3,:)) < pmax(3);
down = min(p(1,:)) > pmin(1) && min(p(2,:)) > pmin(2) && min(p(3,:)) > pmin(3); 

inbounds = up && down;
end