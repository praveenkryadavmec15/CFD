% Semi Explicit method based 2D unsteady NS solver on a uniform staggered grid

clc;
clear;

%% Thermophysical Properties %%
Re  = 100;
rho = 1;
mu  = 1/Re;

%% Geometrial Parameters %%

L1 = 1; % Length of the cavity
L2 = 1; % Height of the cavity

imax = input('Enter Grid Size (imax): ');
jmax = imax;

Dx = L1/(imax-2); % Width of the CV
Dy = L2/(jmax-2); % Height of the CV

x = 0:Dx:L1;
y = 0:Dy:L2;

xc(2:imax-1) = (x(2:imax-1) + x(1:imax-2))/2;
xc(1) = x(1); xc(imax) = x(imax-1);
dx(1:imax-1) = xc(2:imax) - xc(1:imax-1);

yc(2:jmax-1) = (y(2:jmax-1) + y(1:jmax-2))/2;
yc(1) = y(1); yc(jmax) = y(jmax-1);
dy(1:jmax-1) = yc(2:jmax) - yc(1:jmax-1);

%% Initial and Boundary Condition %%

u(imax-1,jmax) = 0;
v(imax,jmax-1) = 0;
P(imax,jmax)   = 0;
P_c(imax,jmax) = 0;

Uo = 1; % Lid velocity
u(2:imax-2,jmax) = Uo;  % Top Boundary

%% Governing Parameters %%
est = 1e-3; % Steady state convergence tolerance
e = 1e-8; % tolerance for mass conservation practically equal to zero


%% Predictor Step %%

n = 0;
Error1 = 1;

while Error1 >= est

    % Timestep %

    Dt_adv = (1/(max(max(abs(u)))/Dx + max(max(abs(v)))/Dy));

    Dt_diff = 0.5/((mu/rho)*((1/Dx)^2 + (1/Dy)^2));

    Dt = 0.4 * min(Dt_diff,Dt_adv);

    n = n + 1;

    u_p   = u;  % u predicted
    v_p   = v;  % v predicted
    u_old = u;  % u previous timestep
    v_old = v;
    P_old = P;
    P_p   = P;  % P predicted

    %Non-Dirichlet BCs for Pressure 
    P(:,jmax) = P(:,jmax-1);     % Top Boundary
    P(:,1) = P(:,2);             % Bottom Boundary
    P(1,:) = P(2,:);             % Left Boundary
    P(imax,:) = P(imax-1,:);     % Right Boundary

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
            u_p(i,j) = u_old(i,j) + ((Dt/(rho*Dx*Dy))*(D_u(i,j) - A_u(i,j) + S_u(i,j)));
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
            v_p(i,j) = v_old(i,j) + Dt/(rho*Dx*Dy) * (D_v(i,j) - A_v(i,j) + S_v(i,j));
            m_y_p(i,j) = rho*v_p(i,j);            
        end
    end   

    % Predicted mass flux on the Boundaries
    m_x_p(1,2:jmax-1) = rho*u(1,2:jmax-1); 
    m_x_p(imax-1,2:jmax-1) = rho*u(imax-1,2:jmax-1);

    m_y_p(2:imax-1,1) = rho*v(2:imax-1,1); 
    m_y_p(2:imax-1,jmax-1) = rho*v(2:imax-1,jmax-1);

    %% Mass Source prediction %%
    
    for i = (2:imax-1)
        for j = (2:jmax-1)
            S_m_p(i,j) = (m_x_p(i,j) - m_x_p(i-1,j))*Dy  + (m_y_p(i,j) - m_y_p(i,j-1))*Dx;
        end
    end

    if max(max(S_m_p)) > e
         
    %% Corrector step %%

    P_c = zeros(imax,jmax);
    N = 0;
    while true
        N = N + 1;
                
        P_p_old   = P_p;
        u_p_old   = u_p;
        v_p_old   = v_p;
        m_x_p_old = m_x_p;
        m_y_p_old = m_y_p;
        S_m_p_old = S_m_p;
        
        P_c(:,jmax) = P_c(:,jmax-1); % Top Boundary
        P_c(:,1)    = P_c(:,2);      % Bottom Boundary
        P_c(1,:)    = P_c(2,:);      % Left Boundary
        P_c(imax,:) = P_c(imax-1,:); % Right Boundary

        P_c_old = P_c;
        
        for i = (2:imax-1)
            for j = (2:jmax-1)
                
                A = -Dt*Dy*((P_c_old(i+1,j) - P_c_old(i,j))/dx(i) - ((P_c_old(i,j) - P_c(i-1,j))/dx(i-1))) ;
                S_m_c = A -Dt*Dx*((P_c_old(i,j+1) - P_c_old(i,j))/dy(j) - ((P_c_old(i,j) - P_c(i,j-1))/dy(j-1)));

                aP = Dt*((1/dx(i) + 1/dx(i-1))*Dy + (1/dy(j) + 1/dy(j-1))*Dx);

                P_c(i,j) = P_c_old(i,j) - ((1/aP) * (S_m_p_old(i,j) + S_m_c));

             end
        end

        for i = (1:imax-1)
            for j = (2:jmax-1)
    
                m_x_c(i,j) = -(Dt/dx(i))*(P_c(i+1,j) - P_c(i,j));
                m_x_p(i,j) = m_x_p_old(i,j) + m_x_c(i,j);
                u_p(i,j) = (m_x_p(i,j))/rho;
            end
        end

        for i = (2:imax-1)
            for j = (1:jmax-1)
    
                m_y_c(i,j) = -(Dt/dy(j))*(P_c(i,j+1) - P_c(i,j));
                m_y_p(i,j) = m_y_p_old(i,j) + m_y_c(i,j);
                v_p(i,j) = (m_y_p(i,j))/rho;
            end
        end

        for i = (2:imax-1)
            for j = (2:jmax-1)
                
                S_m_p(i,j)=(m_x_p(i,j) - m_x_p(i-1,j))*Dy + (m_y_p(i,j) - m_y_p(i,j-1))*Dx;
                P_p(i,j) = P_p_old(i,j) + P_c(i,j);
            end
        end
        %disp(max(max(S_m_p)))
        Error2 = max(max((S_m_p)));
        if Error2 <= e
            break
        end
    end
    end
   
    P  = P_p;
    u  = u_p;
    v  = v_p;
    m_x = m_x_p;
    m_y = m_y_p;

    U =u./Uo;
    U_old=u_old./Uo;
    V = v./Uo;
    V_old = v_old./Uo;
    dt=(Uo*Dt)/L2;
    Error1 = max(max(max(abs(U-U_old)/dt)),max(max(abs(V-V_old)/dt)));
    %disp(Error1)
end

figure(1)
[X,Y]=meshgrid(x,yc);
contourf(X,Y,u',15)
colormap turbo
colorbar
title('U-Velocity Field','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')

figure(10)
[X,Y]=meshgrid(xc,y);
contourf(X,Y,v',15)
colormap turbo
colorbar
title('V-Velocity Field','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')

figure(2)
quiver(x,y,u(1:imax-1,1:jmax-1)',v(1:imax-1,1:jmax-1)',"b", "filled","LineWidth",1)
xlim([0 1])
ylim([0 1])
title('Reference Vector: U = 1','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')
exportgraphics(figure(2),"Vel_vector.tif","Resolution",500)

figure(3)
[X,Y]=meshgrid(xc,yc);
contourf(xc, yc, P',15)
colormap turbo
colorbar
title('Pressure Field','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')

figure(4)
y_com = [0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1];
u_com = [0,-0.03717,-0.04192,-0.04775,-0.063434,-0.1015,-0.15662,-0.2109,-0.20581,-0.136641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1];
scatter(u_com,y_com,'filled','s','DisplayName','Ghia et al. (1982)')
hold on
if rem(imax,2) == 0
    Uc = (u(imax/2,:) + u((imax/2)+1,:))/2;
    plot( Uc,yc,'-','DisplayName', sprintf('%d x %d', imax, jmax))
else
    Uc = u((imax+1)/2,:);
    plot(Uc,yc,'-','DisplayName', sprintf('%d x %d', imax, jmax))
end
legend('-DynamicLegend', "FontSize", 14);
title('Centreline Velocity','FontSize',20,'FontWeight','bold')
xlabel('U','FontSize',18,'FontWeight','bold')
ylabel('Y','FontSize',18,'FontWeight','bold')

figure(5)
x_com=[0;0.0625;0.0703;0.0781;0.0938;0.1563;0.2266;0.2344;0.5;0.8047;0.8594;0.9063;0.9453;0.9531;0.9609;0.9688;1];
v_com=[0;0.09233;0.10091;0.1089;0.12317;0.16077;0.17507;0.17527;0.05454;-0.24533;-0.22445;-0.16914;-0.10313;-0.08864;-0.07391;-0.05906;0];
scatter(x_com,v_com,'filled','v','DisplayName','Ghia et al. (1982)')
hold on
if rem(jmax,2) == 0
    Vc = (v(:,jmax/2) + v(:,(jmax/2)+1))/2;
    plot(xc,Vc,'-','DisplayName', sprintf('%d x %d', imax, jmax))
else
    Vc = v(:,(jmax+1)/2);
    plot(xc,Vc,'-','DisplayName', sprintf('%d x %d', imax, jmax))
end
legend('-DynamicLegend', "FontSize", 14);
title('Centreline Velocity','FontSize',20,'FontWeight','bold')
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('V','FontSize',18,'FontWeight','bold')



