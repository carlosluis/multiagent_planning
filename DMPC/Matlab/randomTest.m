function [po,pf] = randomTest(N,pmin,pmax,rmin)

%Generate initial points
po(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
pass = false;

for n = 2:N
    while(~pass)
        candidate = (pmin + (pmax-pmin).*rand(1,3))';
        diff = po - candidate;
        dist = sqrt(sum(diff.^2,1));

        if(dist > rmin)
            po(:,n) = candidate;
            pass = true;   
        end
    end
    pass = false;
end
po = reshape(po,1,3,N);
%Generate final points

pf(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
pass = false;

for n = 2:N
    while(~pass)
        candidate = (pmin + (pmax-pmin).*rand(1,3))';
        diff = pf - candidate;
        dist = sqrt(sum(diff.^2,1));

        if(dist > rmin)
            pf(:,n) = candidate;
            pass = true;   
        end
    end
    pass = false;
end
pf = reshape(pf,1,3,N);