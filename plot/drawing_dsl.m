clc
close all
clear all

% Initial locations on a 2x2 grid

po1 = [-1.0; 1.0];
po2 = [ -0.5; 1.0];
po3 = [ 0.0; 1.0];
po4 = [ 0.5; 1.0];
po5 = [ 1.0; 1.0];
po6 = [ -1.0; 0.5];
po7 = [ -0.5; 0.5];
po8 = [ 0.0; 0.5];
po9 = [ 0.5; 0.5];
po10 = [ 1.0; 0.5];
po11 = [ -1.0; 0.0];
po12 = [ -0.5; 0.0];
po13 = [ 0.0; 0.0];
po14 = [ 0.5; 0.0];
po15 = [ 1.0; 0.0];
po16 = [ -1.0; -0.5];
po17 = [ -0.5; -0.5];
po18 = [ 0.0; -0.5];
po19 = [ 0.5; -0.5];
po20 = [ 1.0; -0.5];
po21 = [-1.0; -1.0];
po22 = [-0.5; -1.0];
po23 = [0.0; -1.0];
po24 = [0.5; -1.0];
po25 = [1.0; -1.0];

po = cat(2,po1,po2,po3,po4,po5,po6,po7,po8,po9,...
    po10,po11,po12,po13,po14,po15,po16,po17,po18,...
    po19,po20,po21,po22,po23,po24,po25);

% D
pf25 = [-1.25;0.5];
pf15 = [-1.25;0.25];
pf22 = [-1.25;0.0];
pf10 = [-1.25;-0.25];
pf5 = [-1.25;-0.5];
pf21 = [-1.0;0.4];
pf4 = [-1.0;-0.4];
pf23 = [-0.75;0.15];
pf24 = [-0.75;-0.15];

% S
pf18 = [0.25;0.5];
pf9 = [-0.25;-0.5];
pf13 = [0.0;0.0];
pf19 = [0.0;0.45];
pf20 = [-0.25;0.35];
pf14 = [-0.25;0.1];
pf8 = [0.0;-0.45];
pf12 = [0.25;-0.1];
pf7 = [0.25;-0.35];

% L
pf16 = [0.75;0.5];
pf17 = [0.75;0.25];
pf11 = [0.75;0.0];
pf6 = [0.75;-0.25];
pf3 = [0.75;-0.5];
pf2 = [1.0;-0.5];
pf1 = [1.25;-0.5];

pf = cat(2,pf1,pf2,pf3,pf4,pf5,pf6,pf7,pf8,pf9,...
    pf10,pf11,pf12,pf13,pf14,pf15,pf16,pf17,pf18,...
    pf19,pf20,pf21,pf22,pf23,pf24,pf25);
N = size(po,2);
colors = distinguishable_colors(N);

figure(1)

for i=1:N
    if i==13
       plot(po(1,i),po(2,i),'pk','LineWidth',1,...
           'MarkerFaceColor',colors(i,:),'markers',15);
    else
       plot(po(1,i),po(2,i),'ok','LineWidth',1,...
           'MarkerFaceColor',colors(i,:),'markers',15);
    end
   text(po(1,i)-0.05,po(2,i)+0.2,num2str(i),'Fontsize',14)        
   set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
   xticks([-1.5  1.5]);
   yticks([-1.5  1.5]);
   set(gca,'FontSize',20)
   xlabel('x [m]')
   ylabel('y [m]')
   box on
   hold on;
   xlim([-1.7,1.7]);
   ylim([-1.7,1.7]);
end
xh = get(gca,'xlabel'); % handle to the label object
p = get(xh,'position'); % get the current position property
p(2) = p(2)/1.2 ;        % double the distance, 
                       % negative values put the label below the axis
set(xh,'position',p)   % set the new position
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = p(1)/1.2 ;        % double the distance, 
                       % negative values put the label below the axis
set(yh,'position',p)   % set the new position

rectangle('Position',[-1.5,-1.5,3,3],'LineStyle','--',...
    'LineWidth',4,'EdgeColor','r')
set(gcf,'color','w');
%%
figure(2)
for i=1:N
   if i==13
       plot(pf(1,i),pf(2,i),'pk','LineWidth',1,...
           'MarkerFaceColor',colors(i,:),'markers',15);
    else
       plot(pf(1,i),pf(2,i),'dk','LineWidth',1,...
           'MarkerFaceColor',colors(i,:),'markers',15);
    end
            
   set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
   xticks([-1.5  1.5]);
   yticks([-1.5  1.5]);
   set(gca,'FontSize',20)
   xlabel('x [m]')
   ylabel('y [m]')
   box on
   hold on;
   xlim([-1.7,1.7]);
   ylim([-1.7,1.7]);
end

% Text labels for DSL are not straightforward
text(pf(1,1)-0.05,pf(2,1)-0.25,num2str(1),'Fontsize',14);
text(pf(1,2)-0.05,pf(2,2)-0.25,num2str(2),'Fontsize',14);
text(pf(1,3)-0.05,pf(2,3)-0.25,num2str(3),'Fontsize',14);
text(pf(1,4)-0.05,pf(2,4)-0.25,num2str(4),'Fontsize',14);
text(pf(1,5)-0.05,pf(2,5)-0.25,num2str(5),'Fontsize',14);
text(pf(1,6)+0.1,pf(2,6),num2str(6),'Fontsize',14);
text(pf(1,7)+0.1,pf(2,7),num2str(7),'Fontsize',14);
text(pf(1,8)-0.05,pf(2,8)-0.25,num2str(8),'Fontsize',14);
text(pf(1,9)-0.05,pf(2,9)-0.25,num2str(9),'Fontsize',14);
text(pf(1,10)-0.25,pf(2,10)-0.1,num2str(10),'Fontsize',14);
text(pf(1,11)+0.1,pf(2,11),num2str(11),'Fontsize',14);
text(pf(1,12)+0.1,pf(2,12),num2str(12),'Fontsize',14);
text(pf(1,13)-0.1,pf(2,13)+0.2,num2str(13),'Fontsize',14);
text(pf(1,14)-0.1,pf(2,14)-0.25,num2str(14),'Fontsize',14);
text(pf(1,15)-0.25,pf(2,15)-0.1,num2str(15),'Fontsize',14);
text(pf(1,16)+0.1,pf(2,16),num2str(16),'Fontsize',14);
text(pf(1,17)+0.1,pf(2,17),num2str(17),'Fontsize',14);
text(pf(1,18)-0.1,pf(2,18)+0.2,num2str(18),'Fontsize',14);
text(pf(1,19)-0.1,pf(2,19)+0.2,num2str(19),'Fontsize',14);
text(pf(1,20)-0.1,pf(2,20)+0.2,num2str(20),'Fontsize',14);
text(pf(1,21)-0.1,pf(2,21)+0.2,num2str(21),'Fontsize',14);
text(pf(1,22)-0.25,pf(2,22)-0.1,num2str(22),'Fontsize',14);
text(pf(1,23)-0.1,pf(2,23)+0.2,num2str(23),'Fontsize',14);
text(pf(1,24)-0.05,pf(2,24)-0.25,num2str(24),'Fontsize',14);
text(pf(1,25)-0.1,pf(2,25)+0.2,num2str(25),'Fontsize',14);


xh = get(gca,'xlabel'); % handle to the label object
p = get(xh,'position'); % get the current position property
p(2) = p(2)/1.2 ;        % double the distance, 
                       % negative values put the label below the axis
set(xh,'position',p)   % set the new position
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = p(1)/1.2 ;        % double the distance, 
                       % negative values put the label below the axis
set(yh,'position',p)   % set the new position

rectangle('Position',[-1.5,-1.5,3,3],'LineStyle','--',...
    'LineWidth',4,'EdgeColor','r')



set(gcf,'color','w');

% Get the minimum distance between agents
idx = 1;
for i = 1:N
    for j = 1:N
        if (i~=j)
            differ = pf(:,i)-pf(:,j);
            dist(idx) = sqrt(sum(differ.^2));
            idx = idx + 1;
        end
    end
end
