function [po,pf] = randomTest(N,pmin,pmax,rmin)
max_iter = 200000;
%Generate initial points

pass = false;

while(~pass)
    po = [];
    po(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
    for n = 2:N
        tries = 0;
        pass = false;
        while(~pass && tries <= max_iter)
            candidate = (pmin + (pmax-pmin).*rand(1,3))';
            diff = po - candidate;
            dist = sqrt(sum(diff.^2,1));

            if(dist > rmin)
                po(:,n) = candidate;
                pass = true;   
            end
            tries = tries + 1;
        end
        if (tries > max_iter)
            break;
        end
    end
end
po = reshape(po,1,3,N);

%Generate final points

pass = false;

while(~pass)
    pf = [];
    pf(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
    for n = 2:N
        tries = 0;
        pass = false;
        while(~pass && tries <= max_iter)
            candidate = (pmin + (pmax-pmin).*rand(1,3))';
            diff = pf - candidate;
            dist = sqrt(sum(diff.^2,1));

            if(dist > rmin)
                pf(:,n) = candidate;
                pass = true;   
            end
            tries = tries + 1;
        end
        if (tries > max_iter)
            break;
        end
    end
end
pf = reshape(pf,1,3,N);