M = dlmread('trajectories.txt','');
N = M(1,1);
K = size(M,2);
pmin = M(1,2:4);
pmax = M(1,5:7);

po = M(2:4,1:N);
pf = M(5:7,1:N);

all_pos = M(8:end,:);
pk = [];

for i=1:N
    pk(:,:,i) = all_pos(3*(i-1)+1:3*i,:);
end

%%
figure(1)
colors = distinguishable_colors(N);

set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
    for k = 1:K
        for i = 1:N
            hold on;
            grid on;
            xlim([pmin(1),pmax(1)])
            ylim([pmin(2),pmax(2)])
            zlim([0,pmax(3)])
            plot3(pk(1,k,i),pk(2,k,i),pk(3,k,i),'o',...
                'LineWidth',2,'Color',colors(i,:));
            plot3(po(1,i), po(2,i), po(3,i),'^',...
                  'LineWidth',2,'Color',colors(i,:));
            plot3(pf(1,i), pf(2,i), pf(3,i),'x',...
                  'LineWidth',2,'Color',colors(i,:));    
        end
    drawnow
    end
    clf
    pause(0.1)
end
%%
