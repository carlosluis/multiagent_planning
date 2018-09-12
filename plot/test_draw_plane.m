p1 = [-1.3,0,1];
p2 = [2,0.8,1.4];
v = p2-p1;
new_p = p2 - v/norm(v)*0.5;
w = null(v); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(-5:5); % Provide a gridwork (you choose the size)
X = new_p(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = new_p(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = new_p(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z)
xlim([-3,3])
ylim([-3,3])
zlim([-3,3])
hold on;
plot3(p1(1),p1(2),p1(3),'ob');
plot3(p2(1),p2(2),p2(3),'or');
plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)])