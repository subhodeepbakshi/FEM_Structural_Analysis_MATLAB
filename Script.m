clear all;
close all;
clc;
disp('-----------------------------------');
syms eta psi u1x u1y u2x u2y u3x u3y u4x u4y u5x u5y u6x u6y u7x u7y u8x u8y u9x u9y u10x u10y;
assume(eta,'real');
assume(psi,'real');

%   Declaring the Nodes present at each node
disp('Nodes present at each Element :');
Node = [1 2 6 5;
    2 3 7 6;
    3 4 10 7;
    5 6 9 8;
    6 7 10 9]

% Declaring the x-y coordinates of each node
disp('Coordinates of each node : ');
Cor = [0 0;
    3.33 0;
    6.66 0;
    9.99 0;
    0 2.5;
    3.33 2.5;
    6.66 2.5;
    0 5;
    3.33 5;
    9.99 5]

% Displacement Boundary conditions
% BC = [1 1 0;
%     1 2 0;
%     2 2 0;
%     3 2 0;
%     4 2 0;
%     5 1 0;
%     8 1 0];
%
% Force Boundary Conditions
% FC = [8 2 (10*3.33/2);
%     9 2 ((10*3.33/2)+(10*6.66/2));
%     10 2 (10*6.66/2);
%     10 1 4.2;
%     4 1 1.8];

K_Global = zeros(20,20);

% Given Data
disp('Elastic Modulus :');
Em = 1e6

disp('Poisson Ratio');
Pr = 0.3

disp('Finding value of E using the elastic modulus and poisson ratio :');
E = (Em/((1+Pr)*(1-2*Pr)))*[1-Pr Pr 0;Pr 1-Pr 0;0 0 (1-2*Pr)/2]

for i = 1:5

    disp(['Working on Element : ' num2str(i)]);
    disp('-------------------');

    X_jacob = 0;
    Y_jacob = 0;

    for j = 1:4
        syms eta psi;
        assume(eta,'real');
        assume(psi,'real');
        eta ;
        psi ;
        N1 = (1/4)*(1-psi-eta+eta*psi);
        N2 = (1/4)*(1-psi+eta-psi*eta);
        N3 = (1/4)*(1+eta+psi+eta*psi);
        N4 = (1/4)*(1-eta+psi-psi*eta);

        N_local = [N1;N2;N3;N4];

        N = [N1 0 N2 0 N3 0 N4 0;
            0 N1 0 N2 0 N3 0 N4];

        disp(['Using the Coordinates of ' num2str(Node(i,j)) 'th Node' ...
            ' ' num2str(Cor(Node(i,j),1)) ' ' num2str(Cor(Node(i,j),2))]);
        % summing up the values in order to find Jacobian Matrix
        X_jacob = X_jacob + N_local(j)*Cor(Node(i,j),1);
        Y_jacob = Y_jacob + N_local(j)*Cor(Node(i,j),2);
    end;

    disp(' ');
    disp('The Shape Functions are :');
    N1
    N2
    N3
    N4

    % Jacobian Matrix
    rho_x_eta = diff(X_jacob,eta);
    rho_x_psi = diff(X_jacob,psi);
    rho_y_eta = diff(Y_jacob,eta);
    rho_y_psi = diff(Y_jacob,psi);

    disp(' ');
    disp(['Jacobian Matrix for the Element ' num2str(i) ' is :']);
    J = [rho_x_eta rho_y_eta;rho_x_psi rho_y_psi]
    disp(['Inverse Jacobian Matrix for the Element ' num2str(i) ' is :']);
    invJ = inv(J)

    % To Build B Matrix

    rho_N1_rho_x = diff(N1,eta)*invJ(1,1) + diff(N1,psi)*invJ(1,2);
    rho_N1_rho_y = diff(N1,eta)*invJ(2,1) + diff(N1,psi)*invJ(2,2);

    rho_N2_rho_x = diff(N2,eta)*invJ(1,1) + diff(N2,psi)*invJ(1,2);
    rho_N2_rho_y = diff(N2,eta)*invJ(2,1) + diff(N2,psi)*invJ(2,2);

    rho_N3_rho_x = diff(N3,eta)*invJ(1,1) + diff(N3,psi)*invJ(1,2);
    rho_N3_rho_y = diff(N3,eta)*invJ(2,1) + diff(N3,psi)*invJ(2,2);

    rho_N4_rho_x = diff(N4,eta)*invJ(1,1) + diff(N4,psi)*invJ(1,2);
    rho_N4_rho_y = diff(N4,eta)*invJ(2,1) + diff(N4,psi)*invJ(2,2);

    B = [rho_N1_rho_x 0 rho_N2_rho_x 0 rho_N3_rho_x 0 rho_N4_rho_x 0;
        0 rho_N1_rho_y 0 rho_N2_rho_y 0 rho_N3_rho_y 0 rho_N4_rho_y;
        rho_N1_rho_y rho_N1_rho_x rho_N2_rho_y rho_N2_rho_x rho_N3_rho_y rho_N3_rho_x rho_N4_rho_y rho_N4_rho_x];


    if  i/3 == 1
        B_3 = B;
    end;

    % To Form I
    I = B' * E * B * det(J);

    % To perform Gauss Quadrature
    g = sqrt(3)/3;
    disp('The Gauss Quadrature points of the Element are : ');
    eta_psi = [-g -g; g -g; g g; -g g]

    % Variable to store the Local Stiffness Matrix
    I_exact = 0;

    for c = 1:4
        eta = eta_psi(c,1);
        psi = eta_psi(c,2);
        I_exact = I_exact + subs(I);
        I_exact = double(I_exact);
    end;
    disp(' ');
    disp(['Local Stiffness Matrix for the ' num2str(i) 'th Element is ']);
    I_exact

    disp('-----------------------------');

    A = [Node(i,:)];

    for m = 1:8
        B = [2*A(1,1)-1 2*A(1,1) 2*A(1,2)-1 2*A(1,2) 2*A(1,3)-1 2*A(1,3) 2*A(1,4)-1 2*A(1,4)];
        if rem(m,2) == 1
            pos = A(1,(m-1)/2 + 1);
            for l = 1: 8
                K_Global(B(l),2*pos-1) = K_Global(B(l),2*pos-1) + I_exact(l,m);
            end;
        elseif rem(m,2) == 0
            pos = A(1,(m)/2);
            for h = 1:8
                K_Global(B(h),2*pos) = K_Global(B(h),2*pos) + I_exact(h,m);
            end;
        end;
    end;
end;

disp('The Global Stiffness Matrix is : ')
round(vpa(double(K_Global)),3)

disp('The D vector is a 20x1 and is given as :');
d = [u1x;u1y;u2x;u2y;u3x;u3y;u4x;u4y;u5x;u5y;u6x;u6y;u7x;u7y;u8x;u8y;u9x;u9y;u10x;u10y]

disp('The Global Load Vector is a 20x1 and is given as :');
r = [ 0;0;0;0;0;0;1800;0;0;0;0;0;0;0;0;10000*3.33/2;0;10000*9.99/2;4200;10000*6.66/2]

disp('-----------------------------------');


u1x = 0;
u1y = 0;
u2y = 0;
u3y = 0;
u4y = 0;
u5x = 0;
u8x = 0;
displace = subs(d);
displace;

del=[];
q=1;
for h=1:length(displace)
    if displace(h)==0
        del(q)=h;
        q=q+1;
    end;
end;

for z=1:length(del)
    K_Global(del(z),:)=0;
    K_Global(del(z),del(z))=1;
end

KGlobal_reduced=K_Global;
Dsolve=linsolve(KGlobal_reduced,r);

for i=1:length(d)
    d(i)=Dsolve(i);
end

disp('The Nodal displacements (From Node 1 -> Node 10 x,y in order) are - ')
DGlobal = double(round(vpa(d),5));
reshape(DGlobal,[2,10])

% To find the strains and Stresses for Element A ie Element 3
D_B3 = [DGlobal(5) ;DGlobal(6); DGlobal(7); DGlobal(8); DGlobal(19); DGlobal(20); DGlobal(13); DGlobal(14)];
E_strain =  B_3 * D_B3;
eta_psi = [-g -g; g -g; g g; -g g];
disp('To find the strains and Stresses for Element A ie Element 3');
for c = 1:4
    E_strain_eta_psi = 0;
    disp(['Strain and Stress at Gaussian Integration Points eta = ' num2str(eta_psi(c,1)) ' and psi = ' num2str(eta_psi(c,2))]);
    eta = eta_psi(c,1);
    psi = eta_psi(c,2);

    format short g;
    E_strain_eta_psi = double(round(vpa(subs(E_strain)),4))
    Sigma = double(round(vpa(E * E_strain_eta_psi),4))
end;


% Plotting Original and Deformed Structures in the Same Plot
figure;

% Coordinates of each original node
x_nodes = Cor(:, 1);
y_nodes = Cor(:, 2);

% Plotting original nodes
scatter(x_nodes, y_nodes, 'bo', 'filled');
hold on;

% Coordinates of each deformed node
u_x = DGlobal(1:2:end); 
u_y = DGlobal(2:2:end); 

% Plotting deformed nodes
scatter(x_nodes + u_x, y_nodes + u_y, 'ro', 'filled');

% Connecting original nodes with lines
for i = 1:size(Node, 1)
    plot(Cor(Node(i, [1:end 1]), 1), Cor(Node(i, [1:end 1]), 2), 'k-');
end

% Connecting deformed nodes with lines
for i = 1:size(Node, 1)
    plot(Cor(Node(i, [1:end 1]), 1) + u_x(Node(i, [1:end 1])), Cor(Node(i, [1:end 1]), 2) + u_y(Node(i, [1:end 1])), 'b--');
end

title('Comparison between Original and Deformed Structure');
xlabel('X-axis');
ylabel('Y-axis');
legend('Original Node Placement', 'Deformed Node Placement');
grid on;
hold off;