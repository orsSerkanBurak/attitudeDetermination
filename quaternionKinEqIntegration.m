clc;clear;close all;
initialState = [0.408248, 0, 0.408248, 0.816497];
tEnd = 42;
step = 0.01;
[t,sol] = ode45(@quaternionDynamics, 0:step:tEnd, initialState);
norm=sqrt((sol(end,2))^2+(sol(end,3))^2+(sol(end,4))^2)
figure(1)
plot(t,sol,'LineWidth',2);
legend('$\beta_{0}$','$\beta_{1}$','$\beta_{2}$','$\beta_{3}$','Interpreter','latex')
grid on;

function  quaternion = quaternionDynamics(dt,quaternion)
omegaDot = 20*pi/180.*[sin(0.1*dt);0.01;cos(0.1*dt)];
omega = [0;omegaDot];
%b0=quaternion(1);b1=quaternion(2);b2=quaternion(3);b3=quaternion(4);
b0=quaternion(1);b1=quaternion(2);b2=quaternion(3);b3=quaternion(4);
B =0.5*[b0 -b1 -b2 -b3;
        b1 b0 -b3 b2;
        b2 b3 b0 -b1;
        b3 -b2 b1 b0;];
quaternion = B*omega;    
end