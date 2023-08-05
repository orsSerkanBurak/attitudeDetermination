clc;clear;close all;
initialState = pi/180.*[40 , 30 , 80];
tEnd = 42;
step = 0.01;
[t,sol] = ode45(@omegaDot_B, 0:step:tEnd, initialState);
norm=sqrt((sol(end,1))^2+(sol(end,2))^2+(sol(end,3))^2)
figure(1)
plot(t,sol,'LineWidth',2);
legend('$\psi$','$\theta$','$\phi$','Interpreter','latex')
grid on;
function eulerRate = omegaDot_B(dt,~)
eulerRate = 20*pi/180.*[sin(0.1*dt);0.01;cos(0.1*dt)];        
end