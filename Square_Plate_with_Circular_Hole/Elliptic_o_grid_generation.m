% Curvilinear Grid Generation using elliptic grid generation %

L = 10; % Edge length of the square plate 
imax = 10;
jmax = 5;
zetalines = imax-1;
etalines = jmax-1;

e = 1e-4; % Convergence tolerance for grid convergence

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
% Calculating Geometrical Parameters 

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
