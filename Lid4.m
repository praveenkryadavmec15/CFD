clc;
clear;
close all;
%% Thermo physical properties
rho=1;
re=100;
mu=1/re;
nu=mu/rho;
%% Geometric Parameters
l1_L=0;
l1_R=1;
l2=1;
l1=l2;
%% Caluclation of Geometric Parameters
for n=[ 42]
imax=n;
jmax=n;

delx=l1/(imax-2);
dely=l2/(jmax-2);
x(1)=l1_L;
x(imax)=l1_R;
x(2:imax-1)=(l1_L+delx/2:delx:l1_R-delx/2);
dx(1)=delx/2;
dx(imax-1)=dx(1);
dx(2:imax-2)=delx;
y=x;
dy=dx;
end
%% Initial Condition and Boundary Conditions
u_0=1;
v_0=0;
v(1:imax,1:jmax)=0;
u(1:imax,1:jmax-1)=0;
u(1:imax,jmax)=u_0;
p(1:imax,1:jmax)=0;
%% Convergence Tolerance
epsilon_ss=0.001;
convergence_tol=1;
epsilon_mc=10^(-8);
%% Calculation of time step
delt_cfl=((u_0/delx)+(v_0/dely))^(-1);
delt_gfn=((2*nu)^(-1))*(((1/(delx^2))+(1/(dely^2)))^(-1));
delt=min(delt_gfn,delt_cfl);
delt=0.005;
% mmvbc
%% Predictor Velocity
niter = 0;
while convergence_tol>epsilon_ss

niter = niter + 1;
u_star=u;
v_star=v;
u_old=u;
v_old=v;
p_old=p;
p_dash=p;
p_star=p;
%% Prediction of u,v velocities and the mass flux
% u_CV calculations
for i=1 : imax-2
for j=2 : jmax-1
mux(i,j)=0.5*rho*(u(i+1,j)+u(i,j));
muxp(i,j)=max(mux(i,j),0);
muxm(i,j)=min(mux(i,j),0);
aux(i,j)=(muxp(i,j)*u(i,j))+(muxm(i,j)*u(i+1,j));
dux(i,j)=mu*((u(i+1,j)-u(i,j))/delx);
end
end
for i=2:imax-2
for j=1:jmax-1
muy(i,j)=0.5*rho*(v(i+1,j)+v(i,j));
muyp(i,j)=max(muy(i,j),0);
muym(i,j)=min(muy(i,j),0);
auy(i,j)=(muyp(i,j)*u(i,j))+(muym(i,j)*u(i,j+1));
if j==1||j==jmax-1
duy(i,j)=2*mu*((u(i,j+1)-u(i,j))/dely);
else
duy(i,j)=mu*((u(i,j+1)-u(i,j))/dely);
end
end
end
% v-CV calculations
for i=1 : imax-1
for j=2 : jmax-2
mvx(i,j)=0.5*rho*(u(i,j+1)+u(i,j));
mvxp(i,j)=max(mvx(i,j),0);
mvxm(i,j)=min(mvx(i,j),0);
avx(i,j)=(mvxp(i,j)*v(i,j))+(mvxm(i,j)*v(i+1,j));
if i==1||j==imax-1
dvx(i,j)=2*mu*((v(i+1,j)-v(i,j))/delx);
else
dvx(i,j)=mu*((v(i+1,j)-v(i,j))/delx);
end
end
end
for i=2:imax-1
for j=1:jmax-2
mvy(i,j)=0.5*rho*(v(i,j+1)+v(i,j));
mvyp(i,j)=max(mvy(i,j),0);
mvym(i,j)=min(mvy(i,j),0);
avy(i,j)=(mvyp(i,j)*v(i,j))+(mvym(i,j)*v(i,j+1));
dvy(i,j)=mu*((v(i,j+1)-v(i,j))/dely);
end
end
for i=2:imax-2
for j=2:jmax-1
Aux(i,j)=((aux(i,j)-aux(i-1,j))*dely)+((auy(i,j)-auy(i,j-1))*delx);
Dux(i,j)=((dux(i,j)-dux(i-1,j))*dely)+((duy(i,j)-duy(i,j-1))*delx);
Su(i,j)=(p(i,j)-p(i+1,j))*dely;
u_star(i,j)=u(i,j)+((delt/(rho*delx*dely))*(Dux(i,j)-Aux(i,j)+Su(i,j)));
mx_star(i,j)=rho*(u_star(i,j));
end
end
for i=2:imax-1
for j=2:jmax-2
Avx(i,j)=((avx(i,j)-avx(i-1,j))*dely)+((avy(i,j)-avy(i,j-1))*delx);
Dvx(i,j)=((dvx(i,j)-dvx(i-1,j))*dely)+((dvy(i,j)-dvy(i,j-1))*delx);
Sv(i,j)=(p(i,j)-p(i,j+1))*delx;
v_star(i,j)=v(i,j)+((delt/(rho*delx*dely))*(Dvx(i,j)-Avx(i,j)+Sv(i,j)));
my_star(i,j)=rho*(v_star(i,j));
end
end
 
for j=1:jmax
mx_star(1,j)=rho*u(1,j);
mx_star(imax-1,j)=rho*u(imax-1,j);
end
for i=1:imax
my_star(i,1)=rho*v(i,1);
my_star(i,jmax-1)=rho*v(i,jmax-1);
end
for i=2:imax-1
for j=2:jmax-1
Sm_star(i,j)=((mx_star(i,j)-mx_star(i-1,j))*dely)+((my_star(i,j)-my_star(i,j-1))*delx);
end
end
p_dash(1:imax,1:jmax)=0;
while max(max(Sm_star))>epsilon_mc
p_dash(1,:)=p_dash(2,:);
p_dash(imax,:)=p_dash(imax-1,:);
p_dash(:,1)=p_dash(:,2);
p_dash(:,jmax)=p_dash(:,jmax-1);
%p_dash(imax,jmax)=0;
p_dash_old=p_dash;
for i=2:imax-1
for j=2:jmax-1
ap=delt*((((1/dx(i))+(1/dx(i-1)))*dely)+(((1/dy(j))+(1/dy(j-1)))*delx));
sdashm=((-delt*dely)*((((p_dash_old(i+1,j)-p_dash_old(i,j))/dx(i))))-(-delt*dely)*((p_dash_old(i,j)-p_dash(i-1,j))/dx(i-1)))+((-delt*delx)*((((p_dash_old(i,j+1)-p_dash_old(i,j))/dy(j)))-((p_dash_old(i,j)-p_dash(i,j-1))/dy(j-1))));
p_dash(i,j)=p_dash_old(i,j)-((Sm_star(i,j)+sdashm)/ap);
end
end
for i=1:imax-1
for j=2:jmax-1
mx_dash(i,j)=-(delt/dx(i))*(p_dash(i+1,j)-p_dash(i,j));
mx_star(i,j)=mx_star(i,j)+mx_dash(i,j);
u_star(i,j)=mx_star(i,j)/rho;
end
end
for i=2:imax-1
for j=1:jmax-1
my_dash(i,j)=-(delt/dy(j))*(p_dash(i,j+1)-p_dash(i,j));
my_star(i,j)=my_star(i,j)+my_dash(i,j);
v_star(i,j)=my_star(i,j)/rho;
end
end
for i=2:imax-1
for j=2:jmax-1
Sm_star(i,j)=((mx_star(i,j)-mx_star(i-1,j))*dely)+((my_star(i,j)-my_star(i,j-1))*delx);
p_star(i,j)=p_star(i,j)+p_dash(i,j);
end
end
end
 
p=p_star;
u=u_star;
v=v_star;
mx=mx_star;
my=my_star;
U=u./u_0;
U_old=u_old./u_0;
V=v./u_0;
V_old=v_old./u_0;
dt=(u_0*delt)/l2;
convergence_tol=max(max(max(abs(U-U_old)/dt)),max(max(abs(V-V_old)/dt)))
% figure on
% ut=transpose(u)
% contourf(ut)
% hold on
end

figure(1)
contourf(u')