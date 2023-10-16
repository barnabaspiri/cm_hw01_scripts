close all
clear
clc

disp('---- PRELIMINARY CALCULATIONS START ----')

% The undeformed and deformed nodal coordinates
XY = [ [-2,0]; [3,1]; [4,3]; [-2,5] ];
xy = [ [6,9]; [7,1]; [10,4]; [8,8] ];

% The displacement vector of the nodes
U = xy - XY;

% Vector of shape functions
syms xi eta
N = 1/4*[ (1-xi)*(1-eta), ...
          (1+xi)*(1-eta), ...
          (1+xi)*(1+eta), ...
          (1-xi)*(1+eta) ];

% Coordinates of point G in undeformed and deformed configurations
G = [ ssubs(dot(XY(:,1), N), xi, 0, eta, 0), ...
           ssubs(dot(XY(:,2), N), xi, 0, eta, 0) ]

g = [ ssubs(dot(xy(:,1), N), xi, 0, eta, 0), ...
           ssubs(dot(xy(:,2), N), xi, 0, eta, 0) ]

% The displacement in xi and eta
u_xi = dot(U(:,1), N);
u_eta = dot(U(:,2), N);

% Interpolation of the reference coordinates
x = dot(XY(:,1), N);
y = dot(XY(:,2), N);

% Seeing which solution pair is the correct one
syms X Y
sol1 = solve([x==X, y==Y], [xi, eta]);

sol11 = [ssubs(sol1.xi(1), X, G(1), Y, G(2)), ...
         ssubs(sol1.eta(1), X, G(1), Y, G(2))]

sol12 = [ssubs(sol1.xi(2), X, G(1), Y, G(2)), ...
         ssubs(sol1.eta(2), X, G(1), Y, G(2))]

% Solution12 is the right mapping, thus the motion is given as x = chi(X)
chi = [ssubs(u_xi, xi, sol1.xi(2), eta, sol1.eta(2)), ...
       ssubs(u_eta, xi, sol1.xi(2), eta, sol1.eta(2))] + [X,Y];

% Check if the mapping of point G is right
G_to_g = ssubs(chi,X,G(1),Y,G(2))

% The inverse motion is given as X = chi_inv(x)
syms x y
sol2 = solve([chi(1)==x, chi(2)==y], [X,Y]);

% Seeing which solution pair is correct: plugging back g and the result should be G
sol21 = eval([ssubs(sol2.X(1), x, g(1), y, g(2)), ...
              ssubs(sol2.Y(1), x, g(1), y, g(2))])

sol22 = eval([ssubs(sol2.X(2), x, g(1), y, g(2)), ...
              ssubs(sol2.Y(2), x, g(1), y, g(2))])

% Solution21 is the right mapping, thus chi_inv
chi_inv = [sol2.X(1), ...
           sol2.Y(1)];

% Check if the mapping of point G is right
g_to_G = ssubs(chi_inv,x,g(1),y,g(2))

disp('---- PRELIMINARY CALCULATIONS END ----')
fprintf('\n\n\n')

%% Task 1. Position of point G in the reference and spatial configurations
disp('---- TASK 1 ----')

G = eval(G)
g = eval(g)

%% Task 2. The displacement vector for point G
disp('---- TASK 2 ----')

U_G = g - G

%% Task 3. Deformation gradient for point G
disp('---- TASK 3 ----')

% Computing the deformation gradient in the material field
F = [[diff(chi(1),X), diff(chi(1),Y), 0]; ...
     [diff(chi(2),X), diff(chi(2),Y), 0]; ...
     [      0,               0,       1]];

% Computing the deformation gradient at point G
F_G = eval(ssubs(F,X,G(1),Y,G(2)))

%% Task 4. Green-Lagrange and Euler-Almansi strain tensors for point G
disp('---- TASK 4 ----')

% Computing the Right Cauchy-Green deformation tensor
C = F_G'*F_G;

% Computing the Green-Lagrange strain tensor in the material field
E = 1/2*(C - eye(3))

% Computing the inverse of the deformation gradient in the spatial field
F_inv = [[diff(chi_inv(1),x), diff(chi_inv(1),y), 0]; ...
         [diff(chi_inv(2),x), diff(chi_inv(2),y), 0]; ...
         [         0,                  0,         1]];

F_inv_G = eval(ssubs(F_inv,x,g(1),y,g(2)));

% Computing the Cauchy deformation tensor
c = F_inv_G'*F_inv_G;

% Computing the Euler-Almansi strain tensor in the spatial field
e = 1/2*(eye(3) - c)


%% Task 5. Volume ratio and volume strain for point G
disp('---- TASK 5 ----')

% Volume ratio
J = det(F_G)

% Volume strain
eps_V = J - 1

%% Task 6. Isochoric and volumetric part of the deformation gradient for point G
disp('---- TASK 6 ----')

% Isochoric part
F_G_iso = J^(-1/3)*F_G

% Volumetric part
F_G_vol = J^(1/3)*eye(3)

%% Task 7. Principal stretches for point G
disp('---- TASK 7 ----')

[vec, mu] = eig(C);
lam = sqrt(mu)

%% Task 8. Engineering strain in the direction defined by the straight line between node  and node 4 in the reference configuration for point G
disp('---- TASK 8 ----')

NODE2 = [XY(2,1); XY(2,2); 0];
NODE4 = [XY(4,1); XY(4,2); 0];

NLine = NODE2 - NODE4;
NLine = NLine / norm(NLine);

NLine_stretch = sqrt(NLine'*(C)*NLine);
NLine_ES = NLine_stretch - 1

%% Task 9. Determine the area change ratio da/dA
disp('---- TASK 9 ----')

% Computing the Piola deformation tensor
B = F_inv_G*F_inv_G';

% Normal vector of the areas
N = [0; 0; 1];

% Area change ratio da/dA
da_per_dA = det(F_G)*sqrt(N'*B*N)

%% Task 10. Angle of shear for material line elements with initial orientations E1 and E2
disp('---- TASK 10 ----')

% Vectors defining the directions
N1 = [1; 0; 0];
N2 = [0; 1; 0];

% Stretch ratios associated with the initial orientation N1 and N2
lam_N1 = sqrt(N1'*(C)*N1);
lam_N2 = sqrt(N2'*(C)*N2);

% Computing the angle of shear
alpha = acos((N1'*(C)*N2)/(lam_N1*lam_N2))
alpha = rad2deg(alpha)

gamma = 90 - alpha

%% Task 11. Amount of rigid body rotation associated with material particle at G
disp('---- TASK 11 ----')

% Eigenvalues and eigenvectors of C
[NC, mu] = eig(C);
mu1 = sqrt(mu(1,1));
mu2 = sqrt(mu(2,2));
mu3 = sqrt(mu(3,3));
NC1 = NC(:,1);
NC2 = NC(:,2);
NC3 = NC(:,3);

% Spectral representation of U
U = mu1*NC1*NC1' + mu2*NC2*NC2' + mu3*NC3*NC3';

% Determining the rigid body rotation R
R = F_G*inv(U)

% Amount of rigid body rotation: 
theta = acos((trace(R)-1)/2)
theta = rad2deg(theta)

%% Task 12. Components of the spatial Hencky's strain tensor
disp('---- TASK 12 ----')

% Computing the right stretch tensor V
V = F_G*R';

% Computing the spatial Hencky strain tensor
h = logm(V)

%% Functions
% Declare function for double substitute
function substitute = ssubs(expr,old1,new1,old2,new2) 
    substitute = subs(subs(expr,old1,new1),old2,new2);
end


