% Semi Explicit method based 2D unsteady NS solver on a uniform staggered grid

clc; clear;

%% Thermophysical Properties %%
Re = 100;
rho = 1;
Cp = 1;
mu = 1/Re;

%% Geometrial Parameters

L1 = 1; % Length of the cavity
L2 = 1; % Height of the cavity
U = 1; % Velocity of the lid

imax = 42;
jmax = 42;

Dx = L1/(imax-2); Dy = L2/(jmax-2);

x = 0:Dx:L1;
y = 0:Dy:L2;

xc(2:imax-1) = (x(2:imax-1) + x(1:imax-2))/2;
xc(1) = x(1); xc(imax) = x(imax-1);
%dx(1:imax-1) = xc(2:imax) - xc(1:imax-1);

yc(2:jmax-1) = (y(2:jmax-1) + y(1:jmax-2))/2;
yc(1) = y(1); yc(jmax) = y(jmax-1);
%dy(1:jmax-1) = yc(2:jmax) - yc(1:jmax-1);

%% Initial Condition
u_i = 0;
v_i = 0;
P_i = 0;
P_c = 0;

%% Boundary Condition %%
U = 1; % Lid velocity


%% Governing Parameters %%
%Re = rho*U*L1/mu; % Reynolds Number
est = 1e-3; % Steady state convergence tolerance
e = 1e-8; % tolerance for mass conservation practically equal to zero

%% Matrix for storing u, v & P values

u = zeros(imax-1,jmax);
v = zeros(imax,jmax-1);
P = zeros(imax,jmax);

u(2:imax-2,jmax) = U;                 % Top Boundary
% v(2:imax-1,jmax-1) = 0;  
% P(2:imax-1,jmax) = P(2:imax-1,jmax-1);
% 
% u(2:imax-2,1) = 0;                    % Bottom Boundary
% v(2:imax-1,1) = 0;
% P(2:imax-1,1) = P(2:imax-1,2);
% 
% u(1,2:jmax-1) = 0;                    % Left Boundary
% v(1,2:jmax-2) = 0;
% P(1,2:jmax-1) = P(2,2:jmax-1);
% 
% u(imax-1,2:jmax-1) = 0;               % Right Boundary
% v(imax,2:jmax-2) = 0;
% P(imax,2:jmax-1) = P(imax-1,2:jmax-1);

u_p = u;
v_p = v;
P_p = P;


% Matrix for mass flux mx at the faces of the P CV

m_x = zeros(imax-1,jmax);
%u = zeros(imax-1,jmax);
m_x(2:imax-2,jmax) = rho*u(2:imax-2,jmax);
%u(2:imax-2,jmax) = U;
m_x(2:imax-2,1) = rho*u(2:imax-2,1);
%u(2:imax-2,1) = 0;
m_x(1,2:jmax-1) = rho*u(1,2:jmax-1);
%u(1,2:jmax-1) = 0; 
m_x(imax-1,2:jmax-1) = rho*u(imax-1,2:jmax-1);
%u(imax-1,2:jmax-1) = 0; 
m_x_p = m_x;


m_y = zeros(imax,jmax-1);
%v = zeros(imax,jmax-1);
m_y(2:imax-1,jmax-1) = rho*v(2:imax-1,jmax-1);
%v(2:imax-1,jmax-1) = 0;
m_y(2:imax-1,1) = rho*v(2:imax-1,1);
%v(2:imax-1,1) = 0;
m_y(1,2:jmax-2) = rho*v(1,2:jmax-2);
%v(1,2:jmax-2) = 0;
m_y(imax,2:jmax-2) = rho*v(imax,2:jmax-2);
%v(imax,2:jmax-2) = 0;
m_y_p = m_y;

n = 0;
Error1 = 1;
while Error1 >= est

    n = n + 1;

    % Non-Dirichlet BCs for Pressure

    P(2:imax-1,jmax) = P(2:imax-1,jmax-1);     % Top Boundary
    P(2:imax-1,1) = P(2:imax-1,2);             % Bottom Boundary
    P(1,2:jmax-1) = P(2,2:jmax-1);             % Left Boundary
    P(imax,2:jmax-1) = P(imax-1,2:jmax-1);     % Right Boundary


    u_old = u;
    v_old = v;
    P_old = P;
    u_p_old = u_p;
    v_p_old = v_p;
    P_p_old = P_p;
    

    %%  Timestep %%

    Dt_adv = (1/(max(max(abs(u)))/Dx + max(max(abs(v)))/Dy));

    Dt_diff = 0.5/((mu/rho)*((1/Dx)^2 + (1/Dy)^2));

    Dt = 0.4 * min(Dt_diff,Dt_adv);

    %% Predictor Step

    %% u-velocity, and mass flux (x-dirn) prediciton %%

    for i = (1:imax-2)
        for j = (2:jmax-1)
 
            m_ux(i,j) = rho*(u_old(i+1,j) + u_old(i,j))/2;
            a_ux(i,j) = max(0,m_ux(i,j))*u_old(i,j) + min(0,m_ux(i,j))*u_old(i+1,j);
            d_ux(i,j) = mu*(u_old(i+1,j) - u_old(i,j))/Dx;

        end
    end

    for i = (2:imax-2)
        for j = (1:jmax-1)

            m_uy(i,j) = rho*(v_old(i+1,j) + v_old(i,j))/2;
            a_uy(i,j) = max(0,m_uy(i,j))*u_old(i,j) + min(0,m_uy(i,j))*u_old(i,j+1);
        
            if j == 1 || j == (jmax-1)
                d_uy(i,j) = mu*(u_old(i,j+1) - u_old(i,j))/(Dy/2);
            else
                d_uy(i,j) = mu*(u_old(i,j+1) - u_old(i,j))/Dy;
            end
        end
    end

    % Calculation of total advection, diffusion and source terms in the x-momentum transport

    for i = (2:imax-2) 
        for j = (2:jmax-1)
            A_u(i,j) = (a_ux(i,j) - a_ux(i-1,j))*Dy + (a_uy(i,j) - a_uy(i,j-1))*Dx;
            D_u(i,j) = (d_ux(i,j) - d_ux(i-1,j))*Dy + (d_uy(i,j) - d_uy(i,j-1))*Dx;
            S_u(i,j) = (P_old(i,j) - P_old(i+1,j))*Dy;
            u_p(i,j) = u_old(i,j) + Dt/(rho*Cp*Dx*Dy)*(D_u(i,j) - A_u(i,j) + S_u(i,j));
            m_x_p(i,j) = rho*u_p(i,j);     % Mass Flux at the centroid of u-velocity CV
        end
    end

    %% v-velocity and mass flux (y-dirn) prediction %%

    for j = (2:jmax-2)
        for i = (1:imax-1)
        
            m_vx(i,j) = rho*(u_old(i,j+1) + u_old(i,j))/2;
            a_vx(i,j) = max(0,m_vx(i,j))*v_old(i,j) + min(0,m_vx(i,j))*v_old(i+1,j);

            if i == 1 || i == (imax-1)
                d_vx(i,j) = mu*(v_old(i+1,j) - v_old(i,j))/(Dx/2);
            else
                d_vx(i,j) = mu*(v_old(i+1,j) - v_old(i,j))/Dx;
            end
        end
    end


    for j = (1:jmax-2)
        for i = (2:imax-1)

            m_vy(i,j) = rho*(v_old(i,j+1) + v_old(i,j))/2;
            a_vy(i,j) = max(0,m_vy(i,j))*v_old(i,j) + min(0,m_vy(i,j))*v_old(i,j+1);
            d_vy(i,j) = mu*(v_old(i,j+1) - v_old(i,j))/Dy;

        end
    end

    % Calculating total advection, diffusion and source term in the y-momentum transport

    for i = (2:imax-1) 
        for j = (2:jmax-2)
    
            A_v(i,j) = (a_vx(i,j) - a_vx(i-1,j))*Dy + (a_vy(i,j) - a_vy(i,j-1))*Dx;
            D_v(i,j) = (d_vx(i,j) - d_vx(i-1,j))*Dy + (d_vy(i,j) - d_vy(i,j-1))*Dx;
            S_v(i,j) = (P_old(i,j) - P_old(i,j+1))*Dx;
            v_p(i,j) = v_old(i,j) + Dt/(rho*Cp*Dx*Dy) * (D_v(i,j) - A_v(i,j) + S_v(i,j));
            m_y_p(i,j) = rho*v_p(i,j);
            
        end
    end    

    %% Mass Source prediction %%
    
    for i = (2:imax-1)
        for j = (2:jmax-1)
            S_m_p(i,j) = (m_x_p(i,j) - m_x_p(i-1,j))*Dy  + (m_y_p(i,j) - m_y_p(i,j-1))*Dx;
        end
    end

    %disp(S_m_p)

    Error2 = max(abs(S_m_p),[],"all");
    
    P_c = zeros(imax,jmax);

    N = 0;
      
    %% Corrector step %%

    while Error2 >= e
    
        N = N + 1;

        %P_c = zeros(imax,jmax); 

        P_c(2:imax-1,jmax) = P_c(2:imax-1,jmax-1); % Top Boundary
        P_c(2:imax-1,1) = P_c(2:imax-1,2);         % Bottom Boundary
        P_c(1,2:jmax-1) = P_c(2,2:jmax-1);         % Left Boundary
        P_c(imax,2:jmax-1) = P_c(imax-1,2:jmax-1); % Right Boundary

        P_c_old = P_c;
  
        for i = (2:imax-1)
            for j = (2:jmax-1)
    
                A = -Dt*((P_c_old(i+1,j) - P_c_old(i,j))/Dx - (P_c_old(i,j) - P_c(i-1,j))/Dx)*Dy;
                S_m_c(i,j) = A -Dt*((P_c_old(i,j+1) - P_c_old(i,j))/Dy - (P_c_old(i,j) - P_c(i,j-1))/Dy)*Dx;
            
                aP = Dt*((1/Dx + 1/Dx)*Dy + (1/Dy + 1/Dy)*Dx);
            
                P_c(i,j) = P_c_old(i,j) - 1/aP * (S_m_p(i,j) + S_m_c(i,j));
            end
        end

        for i = (1:imax-1)
            for j = (2:jmax-1)
    
                m_x_c(i,j) = -Dt*(P_c(i+1,j) - P_c(i,j))/Dx;
                m_x_p(i,j) = m_x_p(i,j) + m_x_c(i,j);
                u_p(i,j) = m_x_p(i,j)/rho;
            end
        end

        for i = (2:imax-1)
            for j = (1:jmax-1)
    
                m_y_c(i,j) = -Dt*(P_c(i,j+1) - P_c(i,j))/Dy;
                %disp(m_y_p(i,j))
                m_y_p(i,j) = m_y_p(i,j) + m_y_c(i,j);
                %disp(m_y_p(i,j))
                v_p(i,j) = m_y_p(i,j)/rho;
            end
        end
    
        for i = (2:imax-1)
            for j = (2:jmax-1)
    
                S_m_p(i,j) = (m_x_p(i,j) - m_x_p(i-1,j))*Dy + (m_y_p(i,j) - m_y_p(i,j-1))*Dx;
                P_p(i,j) = P_p_old(i,j) + P_c(i,j);
            end
        end
        %disp(max(max(S_m_p)))
    
        Error2 = max(max(abs(S_m_p)));

        %m_x_p = m_x_p_new; 
        %m_y_p = m_y_p_new; 
        P_p_old = P_p;
        %S_m_p = S_m_p_new;
    end

    P = P_p;  u = u_p;  v = v_p;  m_x = m_x_p;  m_y = m_y_p;

    Error1 = max(max(max(abs(u(i,j) - u_old(i,j))))/Dt, max(max(abs(v(i,j) - v_old(i,j))))/Dt);

end
%contourf(u)
 figure(1)
 quiver(u(1:imax-1,1:jmax-1)',v(1:imax-1,1:jmax-1)')
% 
% figure(2)
% y_com = [0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1];
% u_com = [0,-0.03717,-0.04192,-0.04775,-0.063434,-0.1015,-0.15662,-0.2109,-0.20581,-0.136641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1];
% %plot(u((imax-2)/2,:),1:jmax,u_3,y_3)
% plot(u_com,y_com, 'k *')
% hold on
% plot(u(imax/2,:),yc)
% legend('Ghia','42x42 grid','22x22 grid','7x7 grid')
% title('Centreline velocity')
% xlabel('U -->')
% ylabel('Y -->')
% 
% 
% x_com=[0;0.0625;0.0703;0.0781;0.0938;0.1563;0.2266;0.2344;0.5;0.8047;0.8594;0.9063;0.9453;0.9531;0.9609;0.9688;1];
% v_com=[0;0.09233;0.10091;0.1089;0.12317;0.16077;0.17507;0.17527;0.05454;-0.24533;-0.22445;-0.16914;-0.10313;-0.08864;-0.07391;-0.05906;0];
% plot(x_com,v_com, 'k *')
% hold on
% plot(v(:,jmax/2),xc, 'g')
% 
% 
% 
% 
% 
% figure(3)
% y=1-(0:Dy:L2);
% y_com = [1;0.9766;0.9688;0.9609;0.9531;0.8516;0.7344;0.6172;0.5;0.4531;0.2813;0.1719;0.1016;0.0703;0.0625;0.0547;0];
% u_com=[1;0.84123;0.78871;0.73722;0.68717;0.23151;0.00332;-0.136641;-0.20581;-0.2109;-0.15662;-0.1015;-0.063434;-0.04775;-0.04192;-0.03717;0];
% plot(u_com,y_com, 'k *')
% xlim([-0.22 1])
% hold on
% plot(ucl,y, 'g')
% ucl31=[1;0.787620129;0.587711723;0.42114596;0.290555896;0.18961776;0.110101319;0.045170674;-0.009748753;-0.057195052;-0.098228157;-0.133016607;-0.161334794;-0.182928957;-0.197740083;-0.206010366;-0.208285895;-0.20535712;-0.198160901;-0.187675667;-0.174825389;-0.160406644;-0.145042783;-0.129160215;-0.112985847;-0.096555039;-0.079727401;-0.062199838;-0.043513792;-0.023049018;0];
% ycl31=[1;0.966666667;0.933333333;0.9;0.866666667;0.833333333;0.8;0.766666667;0.733333333;0.7;0.666666667;0.633333333;0.6;0.566666667;0.533333333;0.5;0.466666667;0.433333333;0.4;
% 0.366666667;0.333333333;0.3;0.266666667;0.233333333;0.2;0.166666667;0.133333333;0.1;0.066666667;0.033333333;0];
% plot(ucl31,ycl31,'r')
% ucl11=[1;0.400693491;0.087599281;-0.073053304;-0.150595989;-0.168701742;-0.150051734;-0.1161176;-0.079522541;-0.042820845;0];
% ycl11=[1;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.1;0];
% plot(ucl11,ycl11,'b')
% legend('Ghia','51x51 grid','31x31 grid','11x11 grid')
% title('Centreline velocity')
% xlabel('U -->')
% ylabel('Y -->')
% figure(3)
% 
% u_plot = rot90(u);
% u_plot = flip(u_plot,1);
% [X,Y]=meshgrid(xc(1:imax-1),yc);
% contourf(X,Y,u_plot,10,"ShowText","on")
% colormap jet;
% colorbar
% title('u-Velocity Field ' ,'FontSize',20,'FontWeight','bold')
% xlabel('X (m)','FontSize',18,'FontWeight','bold') 
% ylabel('Y (m)','FontSize',18,'FontWeight','bold')
% 
% figure(4)
% 
% P_plot = rot90(P);
% P_plot = flip(P_plot,1);
% [X,Y]=meshgrid(xc,yc);
% contourf(X,Y,P_plot,10,"ShowText","on")
% colormap jet;
% colorbar
% title('Pressure Field ' ,'FontSize',20,'FontWeight','bold')
% xlabel('X (m)','FontSize',18,'FontWeight','bold') 
% ylabel('Y (m)','FontSize',18,'FontWeight','bold')
% 
% figure(5)
% 
% [startX,startY] = meshgrid(0.8, 0:0.1:1);
% verts = stream2(xc(2:imax-1), yc(2:jmax-1),u(2:imax-1,2:jmax-1),v(2:imax-1,2:jmax-1),startX,startY);
% lineobj = streamline(verts);
