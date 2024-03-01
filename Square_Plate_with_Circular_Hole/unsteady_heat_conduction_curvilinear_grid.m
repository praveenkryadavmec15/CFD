% Unsteady heat conduction in a square plate with circular hole %

clc; clear;

%% Thermophysical Properties %% 

rho = 7750; % kg/m3 
Cp = 500; % J/(kg.K)
k = 16.2; % W/(m.K)

%% Curvilinear Grid Generation using elliptic grid generation %%

L = 10; % Edge length of the square plate
imax = 10;
jmax = 5;
zetalines = imax-1;
etalines = jmax-1;

e = 1e-4; % Convergence tolerance for grid

R = 2; % Radius of the inner circular hole
theta = linspace(0,-360,imax-1);

Dzeta = 1/(imax-2);
Deta = 1/(jmax-2);

x = ones(imax-1,jmax-1);
y = ones(imax-1,jmax-1);

x(1:imax-1,1) = R*cosd(theta)+L/2; % AB (Inner Boundary, Circular hole)
x(1:imax-1,jmax-1) = [10,10,5,0,0,0,5,10,10]; % DC (Outer Boundary, Square)
x(1,1:jmax-1) = linspace(7,10,jmax-1); % AD (Branch Cut)
x(imax-1,1:jmax-1) = linspace(7,10,jmax-1); % BC (Branch Cut)

y(1:imax-1,1) = R*sind(theta)+L/2; % AB (Inner Boundary, Circular hole) 
y(1:imax-1,jmax-1) = [5,0,0,0,5,10,10,10,5]; % DC (Outer Boundary, Square)
y(1,1:jmax-1) = 5; % AD (Branch Cut)
y(imax-1,1:jmax-1) = 5; % BC (Branch Cut)

N=0;
while true
    N= N + 1;

    x_old = x;
    y_old = y;

    for i = 2:imax-2
        for j = 2:jmax-2

            A(i,j) = (((x_old(i,j+1) - x(i,j-1))/(2*Deta))^2) + (((y_old(i,j+1) - y(i,j-1))/(2*Deta))^2);
            C(i,j) = (((x_old(i+1,j) - x(i-1,j))/(2*Dzeta))^2) + (((y_old(i+1,j) - y(i-1,j))/(2*Dzeta))^2);
            B(i,j) = ((x_old(i+1,j) - x(i-1,j))*(x_old(i,j+1) - x(i,j-1))/(2*Dzeta*2*Deta)) + ((y_old(i+1,j) - y_old(i-1,j))*(y_old(i,j+1) - y_old(i,j-1))/(2*Dzeta*2*Deta));
            
            d_x(i,j) = (x_old(i+1,j+1) + x(i-1,j-1) - x_old(i-1,j+1) - x_old(i+1,j-1))/(4*Deta*Dzeta);
            d_y(i,j) = (y_old(i+1,j+1) + y(i-1,j-1) - y_old(i-1,j+1) - y_old(i+1,j-1))/(4*Deta*Dzeta);
        
            x(i,j) = (A(i,j)*(Deta^2) * (x_old(i+1,j) + x(i-1,j)) + C(i,j)*(Dzeta^2) *(x_old(i,j+1) + x(i,j-1)) - 2*(B(i,j))*(Dzeta^2)*(Deta^2)*(d_x(i,j)))/(2*((A(i,j))*(Deta^2) + (C(i,j))*(Dzeta^2)));
            y(i,j) = (A(i,j)*(Deta^2) * (y_old(i+1,j) + y(i-1,j)) + C(i,j)*(Dzeta^2) *(y_old(i,j+1) + y(i,j-1)) - 2*(B(i,j))*(Dzeta^2)*(Deta^2)*(d_y(i,j)))/(2*((A(i,j))*(Deta^2) + (C(i,j))*(Dzeta^2)));
        end
    end

    error = max(max(max(abs(x - x_old))), max(max(abs(y - y_old))));
    %disp(error)

    if error <= e
        break
    else
        continue
    end
end

figure(1)

plot(x',y','-s',LineWidth=2)
hold on
plot(x,y,'-s',LineWidth=2)

%% Calculating Geometrical Parameters %%

% Centroid of the CV %
for i = 2:imax-1
    for j = 2:jmax-1

        xc(i,j) = (x(i-1,j-1) + x(i-1,j) + x(i,j-1) + x(i,j))/4;
        yc(i,j) = (y(i-1,j-1) + y(i-1,j) + y(i,j-1) + y(i,j))/4;
    end
end

xc(2:imax-1,1) = (x(1:imax-2,1) + x(2:imax-1,1))/2;
xc(2:imax-1,jmax) = (x(1:imax-2,jmax-1) + x(2:imax-1,jmax-1))/2;
xc(1,2:jmax-1) = xc(imax-1,2:jmax-1);
xc(imax,2:jmax-1) = xc(2,2:jmax-1);
xc(1,1) = xc(imax-1,1);   xc(1,jmax) = xc(imax-1,jmax);
xc(imax,1) = xc(2,1);   xc(imax,jmax) = xc(2,jmax);


yc(2:imax-1,1) = (y(1:imax-2,1) + y(2:imax-1,1))/2;
yc(2:imax-1,jmax) = (y(1:imax-2,jmax-1) + y(2:imax-1,jmax-1))/2;
yc(1,2:jmax-1) = yc(imax-1,2:jmax-1);
yc(imax,2:jmax-1) = yc(2,2:jmax-1);
yc(1,1) = yc(imax-1,1);   yc(1,jmax) = yc(imax-1,jmax);
yc(imax,1) = yc(2,1);   yc(imax,jmax) = yc(2,jmax);

scatter(xc',yc','o','filled')
title('Curvilinear Elliptic O Grid','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')
axis square

% Projection of s and t direction along x and y %

for i = 1:imax-1
    for  j = 2:jmax-1

        sx_e_cap(i,j) = (xc(i+1,j) - xc(i,j))/(((xc(i+1,j) - xc(i,j))^2) + ((yc(i+1,j) - yc(i,j))^2))^0.5;
        
        sy_e_cap(i,j) = (yc(i+1,j) - yc(i,j))/(((xc(i+1,j) - xc(i,j))^2) + ((yc(i+1,j) - yc(i,j))^2))^0.5;

        tx_e_cap(i,j) = (x(i,j) - x(i,j-1))/((x(i,j) - x(i,j-1))^2 + (y(i,j) - y(i,j-1))^2)^0.5;
        
        ty_e_cap(i,j) = (y(i,j) - y(i,j-1))/((x(i,j) - x(i,j-1))^2 + (y(i,j) - y(i,j-1))^2)^0.5;

    end
end

for i = 2:imax-1
    for j = 1:jmax-1

        sx_n_cap(i,j) = (xc(i,j+1) - xc(i,j))/(((xc(i,j+1) - xc(i,j))^2) + ((yc(i,j+1) - yc(i,j))^2))^0.5;

        sy_n_cap(i,j) = (yc(i,j+1) - yc(i,j))/(((xc(i,j+1) - xc(i,j))^2) + ((yc(i,j+1) - yc(i,j))^2))^0.5;

        tx_n_cap(i,j) = (x(i,j) - x(i-1,j))/((x(i,j) - x(i-1,j))^2 + (y(i,j) - y(i-1,j))^2)^0.5;

        ty_n_cap(i,j) = (y(i,j) - y(i-1,j))/((x(i,j) - x(i-1,j))^2 + (y(i,j) - y(i-1,j))^2)^0.5;

    end
end

for i = 1:imax-1
    for j = 2:jmax-1

        % Area of the faces of the CV %

        DS_xe(i,j) = y(i,j) - y(i,j-1);
        DS_ye(i,j) = x(i,j-1) - x(i,j);

    end
end

for i = 2:imax-1
    for j = 1:jmax-1

        % Area of the faces of the CV %

        DS_xn(i,j) = -(y(i,j) - y(i-1,j));
        DS_yn(i,j) = x(i,j) - x(i-1,j);

    end
end

DV(1,1:jmax) = 0; DV(imax,1:jmax) = 0; DV(1:imax,1) = 0;  DV(1:imax,jmax) = 0; % Volume of the CV

for i = 1:imax-1
    for j = 2:jmax-1

        dte(i,j) = ((x(i,j) - x(i,j-1))^2 + (y(i,j) - y(i,j-1))^2)^0.5;
    end
end

for i = 2:imax-1
    for j = 1:jmax-1
        dtn(i,j) = ((x(i,j) - x(i-1,j))^2 + (y(i,j) - y(i-1,j))^2)^0.5;
    end
end

for i = 2:imax-1
    for j = 2:jmax-1

        d1(i,j) = (((x(i,j) - x(i-1,j-1))^2) + ((y(i,j) - y(i-1,j-1))^2))^0.5; % Diagonal 1 CV
        d2(i,j) = (((x(i-1,j) - x(i,j-1))^2) + ((y(i-1,j) - y(i,j-1))^2))^0.5; % Diagonal 1 CV

        DV(i,j) = (1/4)*(((2*d1(i,j)*d2(i,j))^2 - (dte(i,j)^2 + dte(i-1,j)^2 - dtn(i,j)^2 - dtn(i,j-1)^2)^2))^0.5; % Vol of CV
    end
end

volume = 0;
for i= 1:imax
    for j = 1:jmax
volume = volume + DV(i,j);
    end
end

for i = 1:imax-1
    for j = 2:jmax-1

        dse(i,j) = ((xc(i+1,j) - xc(i,j))^2 + (yc(i+1,j) - yc(i,j))^2)^0.5;
    end
end

for i = 2:imax-1
    for j = 1:jmax-1

        dsn(i,j) = ((xc(i,j+1) - xc(i,j))^2 + (yc(i,j+1) - yc(i,j))^2)^0.5;
    end
end

for i = 1:imax-1
    for j = 2:jmax-1

        DS_se(i,j) = (DS_xe(i,j)* ty_e_cap(i,j) - DS_ye(i,j) * tx_e_cap(i,j))/(sx_e_cap(i,j)*ty_e_cap(i,j) - sy_e_cap(i,j)*tx_e_cap(i,j));
        
        DS_te(i,j) = (DS_xe(i,j)* sy_e_cap(i,j) - DS_ye(i,j) * sx_e_cap(i,j))/(sx_e_cap(i,j)*ty_e_cap(i,j) - sy_e_cap(i,j)*tx_e_cap(i,j));
        
    end
end

for i = 2:imax-1
    for j = 1:jmax-1

        DS_sn(i,j) = (DS_xn(i,j)* ty_n_cap(i,j) - DS_yn(i,j) * tx_n_cap(i,j))/(sx_n_cap(i,j)*ty_n_cap(i,j) - sy_n_cap(i,j)*tx_n_cap(i,j));
        
        DS_tn(i,j) = (DS_xn(i,j)* sy_n_cap(i,j) - DS_yn(i,j) * sx_n_cap(i,j))/(sx_n_cap(i,j)*ty_n_cap(i,j) - sy_n_cap(i,j)*tx_n_cap(i,j));
    
    end
end

%% Computation of Temperature %%

T = zeros(imax,jmax);
T(1:imax,1) = 1; % AB (Outer Square Boundary)
T(2:imax-1,jmax) = 0; % DC (Outer Square Boundary)
T(1,2:jmax-1) = T(imax-1,2:jmax-1); 
T(imax,2:jmax-1) = T(2,2:jmax-1);

T_old = T;
 
% Governing Parameters %

Dt = 1; % Timestep
Q_gen = 0; % Heat generation
e_ss = 1e-6; % Steady state convergence tolerance for Temperature
    
n = 0;
while true
    n = n + 1;
    t = n*Dt;

    T(1,2:jmax-1) = T(imax-1,2:jmax-1); 
    T(imax,2:jmax-1) = T(2,2:jmax-1);

    
    for i = 1:imax-1
        for j = 1:jmax-1
    
            Tv(i,j) = (DV(i+1,j+1)*T_old(i,j) + DV(i,j)*T_old(i+1,j+1) + DV(i+1,j)*T_old(i,j+1) + DV(i,j+1)*T_old(i+1,j))/(DV(i+1,j+1) + DV(i,j) + DV(i+1,j) + DV(i,j+1));
        end
    end
    
    for i = 1:imax-1
        for j = 2:jmax-1
    
            q_se(i,j) =  -k*(T_old(i+1,j) - T_old(i,j))/dse(i,j);   q_te(i,j) = -k*(Tv(i,j) - Tv(i,j-1))/dte(i,j);
    
        end
    end
    
    for i = 2:imax-1
        for j = 1:jmax-1
    
            q_sn(i,j) = -k*(T_old(i,j+1) - T_old(i,j))/dsn(i,j);    q_tn(i,j) = -k*(Tv(i,j) - Tv(i-1,j))/dtn(i,j);
    
        end
    end
    
    for i = 2:imax-1
        for j = 2:jmax-1
    
            Q_s(i,j) = (q_se(i-1,j)*DS_se(i-1,j) - q_se(i,j)*DS_se(i,j)) + (q_sn(i,j-1)*DS_sn(i,j-1) - q_sn(i,j)*DS_sn(i,j));
            
            Q_t(i,j) = (q_te(i-1,j)*DS_te(i-1,j) - q_te(i,j)*DS_te(i,j)) + (q_tn(i,j-1)*DS_tn(i,j-1) - q_tn(i,j)*DS_tn(i,j));
            
            T(i,j) = T_old(i,j) + Dt/(rho*Cp*DV(i,j))*(Q_s(i,j) + Q_t(i,j) + Q_gen*DV(i,j));
    
        end
    end

    error1 = max(max(abs(T - T_old)))/Dt;
    if rem(n,1000) == 0
      %disp(error1)
    end

    if error1 <= e_ss
        break
    else
        T_old = T;
    end
end

%% Plotting the Temperature Contours %%

figure(2)
contourf(xc',yc',T',10,ShowText="off")
colormap 'turbo'
colorbar
title('Temperature on Cell Centroids','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')
yticks(0:2:10)
axis square

figure(3)
contourf(x',y',Tv',10, ShowText="off")
colormap 'turbo'
colorbar
title('Temperature on Cell Vertices','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')
yticks(0:2:10)
axis square

