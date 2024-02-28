clc
clear
%% Lid Driven cavity flow 
%% Thermophysical properties
L1=1; L2=1; D=1; Re=100; Mu=1/Re; U_o=1;cp=1;
imax=7; jmax=7;
E_ss=10^(-3); E_mc=10^-8;
%% geometrical prameter
dx= L1/(imax-2); dy=L2/(jmax-2);
delx=dx; dely=dy;
%% grid generation
u=zeros(imax-2,jmax);
v=zeros(imax,jmax-2);
p=zeros(imax,jmax);

%% time step
t_adv=1/abs(1/dx); t_diff=(0.5*(1/(dx^2)+1/(dy^2))^-1)*(D/Mu);
dt=0.005*min(t_adv,t_diff);
DT=U_o*dt/L1; %Non dimensional time 
%% Boudary condition

u(2:imax-2,jmax)=1; %Top wall
u(1,2:jmax-1)=0; %left wall
u(imax-1,2:jmax-1)=0; %Rigt wall
u(2:imax-2,1)=0; %Bottom wall 

v(2:imax-1,jmax-1)=0; %Top wall
v(1,2:jmax-2)=0; %Left wall
v(imax,2:jmax-2)=0; %Right wall
v(2:imax-1,1)=0; %Bottom wall 

u_pr=u;
v_pr=v;

p_pr=zeros(imax,jmax);
E_ss_cal=1;
E_mc_cal=1;

%u_pr=u; v_pr=v; p_pr=p;

%% finding fluxes for u control volume-Predictor Step
q=0;
max_iter=100000;
while E_ss_cal>E_ss && q<max_iter
    q=q+1;
    
    for i=1:imax-2
        for j=2:jmax-1
            [m_ux(i,j),a_ux(i,j),d_ux(i,j)]= flux_mx_u(D,dx,Mu,u(i+1,j),u(i,j));
        end
    end
    
    for i=2:imax-2
        for j=1:jmax-1
            [m_uy(i,j),a_uy(i,j),d_uy(i,j)]= flux_my_u(D,Mu,dy,u(i,j+1),u(i,j),v(i+1,j),v(i,j),j,jmax);
        end
    end
    
    for i=2:imax-2
        for j=2:jmax-1
            A_ux(i,j)=(a_ux(i,j)-a_ux(i-1,j))*dy+(a_uy(i,j)-a_uy(i,j-1))*dx;
            D_ux(i,j)=(d_ux(i,j)-d_ux(i-1,j))*dy+(d_uy(i,j)-d_uy(i,j-1))*dx;
            S_u(i,j)=(p(i,j)-p(i+1,j))*dy;
            u_pr(i,j)=u(i,j)+(dt/(D*cp*dx*dy))*(D_ux(i,j)-A_ux(i,j)+S_u(i,j));
            m_x_pr(i,j)=D*u_pr(i,j);
        end
    end
    m_x_pr(1,:)=D*u(1,1:jmax-1);
    m_x_pr(imax-1,:)=D*u(imax-1,1:jmax-1);
    
    %% finding fluxes for v control volume-Predictor Step 
    for i=1:imax-1
        for j=2:jmax-2
            [m_vx(i,j),a_vx(i,j),d_vx(i,j)]=flux_mx_v(D,Mu,dx,u(i,j),u(i,j+1),v(i,j),v(i+1,j),i,imax);
        end
    end
    
    for i=2:imax-1
        for j=1:jmax-2
            [m_vy(i,j),a_vy(i,j),d_vy(i,j)]=flux_my_v(D,Mu,dy,v(i,j+1),v(i,j));
        end
    end
    
    for i=2:imax-1
        for j=2:jmax-2
            A_vx(i,j)=(a_vx(i,j)-a_vx(i-1,j))*dy+(a_vy(i,j)-a_vy(i,j-1))*dx;
            D_vx(i,j)= (d_vx(i,j)-d_vx(i-1,j))*dy+(d_vy(i,j)-d_vy(i,j-1))*dx;
            S_v(i,j)=(p(i,j)-p(i,j+1))*dx;
            v_pr(i,j)= v(i,j)+(dt/(D*dx*dy*cp))*(D_vx(i,j)-A_vx(i,j)+S_v(i,j));
            m_y_pr(i,j)=D*v_pr(i,j);
        end
    end
    
    m_y_pr(:,1)=D*v(1:imax-1,1);
    m_y_pr(:,jmax-1)=D*v(1:imax-1,jmax-1);
    
    %% Mass source predictiom
    
    for i=2:imax-1
        for j=2:jmax-1
            S_m_pr(i,j)=(m_x_pr(i,j)-m_x_pr(i-1,j))*dy+(m_y_pr(i,j)-m_y_pr(i,j-1))*dx;
        end
    end
    
    if S_m_pr<E_mc
        break;
    end
    %% Corrector Step 
    p_cr=zeros(imax,jmax);
    p_pr_new=zeros(imax,jmax);
    max_itt_1=1000;
    w=0;
    while E_mc_cal>E_mc && w<max_itt_1
        w=w+1;
        p_cr(:,jmax)=p_cr(:,jmax-1); %top
        p_cr(1,:)=p_cr(2,:); %Left
        p_cr(:,1)=p_cr(:,2); %Bottom
        p_cr(imax,:)=p_cr(imax-1,:);%Right
        
        for i=2:imax-1
            for j=2:jmax-1          
                ap=dt*(2/delx*dy+2/dely*dx);
                A=-dt*((p_cr(i+1,j)-p_cr(i,j))/delx-((p_cr(i,j)-p_cr(i-1,j))/delx))*dy;
                B=-dt*((p_cr(i,j+1)-p_cr(i,j))/dely-((p_cr(i,j)-p_cr(i,j-1))/dely))*dx;
                S_m_cr=A+B;
                p_cr_new(i,j)=p_cr(i,j)-(1/ap)*(S_m_pr(i,j)-S_m_cr);
            end
        end
        
        p_cr_new(:,jmax)=p_cr_new(:,jmax-1); %top
        p_cr_new(1,:)=p_cr_new(2,:); %Left
        p_cr_new(:,1)=p_cr_new(:,2); %Bottom
        p_cr_new(imax,:)=p_cr_new(imax-1,:);%Right

        u_pr_new(2:imax-2,jmax)=1; %Top wall
        u_pr_new(1,2:jmax-1)=0; %left wall
        u_pr_new(imax-1,2:jmax-1)=0; %Rigt wall
        u_pr_new(2:imax-2,1)=0; %Bottom wall 

        v_pr_new(2:imax-1,jmax-1)=0; %Top wall
        v_pr_new(1,2:jmax-2)=0; %Left wall
        v_pr_new(imax,2:jmax-2)=0; %Right wall
        v_pr_new(2:imax-1,1)=0; %Bottom wall 

   

        for i=1:imax-1
            for j=2:jmax-1
                m_x_cr_new(i,j)=-dt*(p_cr_new(i+1,j)-p_cr_new(i,j))/delx;                
                m_x_pr_new(i,j)=m_x_pr(i,j)+m_x_cr_new(i,j);                
                u_pr_new(i,j)=m_x_pr_new(i,j)/D;               
            end
        end


        for i=2:imax-1
            for j=1:jmax-1                
                m_y_cr_new(i,j)=-dt*(p_cr_new(i,j+1)-p_cr_new(i,j))/dely;               
                m_y_pr_new(i,j)=m_y_pr(i,j)+m_y_cr_new(i,j);                
                v_pr_new(i,j)=m_y_pr_new(i,j)/D;
            end
        end

       w;
       p_pr_new(:,jmax)=p_pr_new(:,jmax-1); %Top             
       p_pr_new(1,:)=p_pr_new(2,:); %left       
       p_pr_new(:,1)=p_pr_new(:,2);%Bottom
       p_pr_new(imax,:)=p_pr_new(imax-1,:); %Right
        
       for i=2:imax-1
           for j=2:jmax-1
               S_m_pr_new(i,j)=(m_x_pr_new(i,j)-m_x_pr_new(i-1,j))*dy+(m_y_pr_new(i,j)-m_y_pr_new(i,j-1))*dx;
               p_pr_new(i,j)=p_pr(i,j)+p_cr_new(i,j);
           end
       end
       E_mc_cal=max(max(abs(S_m_pr_new)));
       % p_pr_new(:,jmax)=p_pr_new(:,jmax-1); %Top             
       % p_pr_new(1,:)=p_pr_new(2,:);         %left       
       % p_pr_new(:,1)=p_pr_new(:,2);         %Bottom
       % p_pr_new(imax,:)=p_pr_new(imax-1,:); %Right
       u_pr=u_pr_new;
       v_pr=v_pr_new;
       S_m_pr=S_m_pr_new;
       m_x_pr=m_x_pr_new;
       m_y_pr=m_y_pr_new;
    end

    DU=max(max(abs(u_pr_new-u)))/DT;
    DV=max(max(abs(v_pr_new-v)))/DT;
    E_ss_cal=max(DU,DV);
    
    u=u_pr_new;
    v=v_pr_new;
    p=p_pr_new;    
    
end

x = 0:L1/(imax-1):L1;
y = 0:L2/(jmax-1):L2;

figure (1)
quiver(u(1:imax-1,1:jmax-1)',v(1:imax-1,1:jmax-1)',x,y)
xlabel('x');
ylabel('y');
title('Velocity vector contours')

figure (2)

y_3 = [0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1];
u_3 = [0,-0.03717,-0.04192,-0.04775,-0.063434,-0.1015,-0.15662,-0.2109,-0.20581,-0.136641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1];

plot (u((imax-2)/2,:),0:L1/(jmax-1):L1,u_3,y_3,'*')
title('U-velocity along verticle Centre line')
legend('42*42','22*22','12*12','7*7','Ghia Re100')
xlabel('U Vel (m/s)');
ylabel('Y');
hold on

figure (3)

x_4=[0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1];
v_4=[0,0.09233,0.10091,0.1089,0.12317,0.16077,0.17507,0.17527,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864,-0.07391,-0.05906,0];
plot (0:L2/(imax-1):L2,v(:,(jmax-2)/2),x_4,v_4,'*')
title('V-velocity along horizontal Centre line')
legend('42*42','22*22','12*12','7*7','Ghia Re100')
xlabel('X');
ylabel('V Vel (m/s)');
hold on
colobar on

figure (4)
contourf(x(1:imax-1),y,u',10)
title('U-velocity contour')
colobar on

figure (5)
Xx = 0:L1/(imax-2):L1;
Yy = 0:L2/(jmax-2):L2;
[X,Y] = meshgrid(Xx,Yy); 
N = 25; 
xstart = max(x)*rand(N,1); 
ystart = max(y)*rand(N,1);
contourf(x,y,p',10)
hold on
h=streamline(X,Y,u(1:41,1:41)',v(1:41,1:41)',xstart,ystart);
set(h,'color','Cyan')
hold off
title('Pressure Distribution & Streamiles')
colobar on
         

%%  THE END OF  :D %%

%% Defined Functions 

function [m,a,d]= flux_mx_u(D,dx,Mu,u_e,u_p)
    m=D*(u_e+u_p)/2;
    a=max(m,0)*u_p+min(m,0)*u_e;
    d=Mu*(u_e-u_p)/dx;
end

function [m,a,d]= flux_my_u(D,Mu,dy,u_n,u_p,v_e,v_p,j,jmax)
    m=D*(v_e+v_p)/2;
    a=max(m,0)*u_p+min(m,0)*u_n;
    if j==1 || j==jmax-1
        d=Mu*(u_n-u_p)/(dy/2);
    else
        d=Mu*(u_n-u_p)/dy;
    end

end

function [m,a,d]=flux_mx_v(D,Mu,dx,u_p,u_n,v_p,v_e,i,imax)
        m=D*(u_n+u_p)/2;
        a=max(m,0)*v_p+min(m,0)*v_e;
        if i==1 || i==imax-1
            d=Mu*(v_e-v_p)/(dx/2);
        else
            d=Mu*(v_e-v_p)/dx;
        end
end

function [m,a,d]=flux_my_v(D,Mu,dy,v_n,v_p)
        m=D*(v_n+v_p)/2;
        a=max(m,0)*v_p+min(m,0)*v_n;
        d=Mu*(v_n-v_p)/dy;        
end






