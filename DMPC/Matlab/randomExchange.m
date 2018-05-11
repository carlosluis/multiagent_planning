function [po,pf] = randomExchange(N,pmin,pmax,rmin)

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

%Make a random permutation of initial states to assign final states

idx_vector = randperm(N);
for n = 1:N
    pf(:,n) = po(:,idx_vector(n));
end

po = reshape(po,1,3,N);
pf = reshape(pf,1,3,N);