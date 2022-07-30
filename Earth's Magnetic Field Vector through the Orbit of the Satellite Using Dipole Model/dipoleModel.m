clc;clear;close all;

n=20;%sequence number

Me = 7.943e15; %the magnetic dipole moment of the Earth in Wb.m
mu = 3.98601e14; %the Earth Gravitational constant in m^3/s^2
i = (80 + (0.5*n))*(pi/180); %inclination converted from deg to rad
we = 7.29e-5; %the spin rate of the Earth in rad/s
epsilon = 11.7*(pi/180); %the magnetic dipole tilt converted to rad
r0 = (6378.14 + 500 + (2*n))*1000; %the distance between the center of mass of the satellite and the Earth in m
w0 = sqrt(mu/r0^3); %the angular velocity of the orbit with respect to the inertial frame in rad/s

dt = 0.1; %the sample time(s)
N = 54000*dt; %iteration

for k = 1:N
%Components of Earth Magnetic Field (Wb/m^2(Tesla))
hx = (Me/r0^3)*((cos(w0*k)*((cos(epsilon)*sin(i))-(sin(epsilon)*cos(i)*cos(we*k)))-(sin(w0*k)*sin(epsilon)*sin(we*k))));
hy = (-(Me/r0^3))*((cos(epsilon)*cos(i))+(sin(epsilon)*sin(i)*cos(we*k)));
hz = ((2*Me)/r0^3)*((sin(w0*k)*((cos(epsilon)*sin(i))-(sin(epsilon)*cos(i)*cos(we*k)))-(2*(sin(w0*k)*sin(epsilon)*sin(we*k)))));

%Components of Direction Cosine Matrix
hx0 = hx/(sqrt(hx^2+hy^2+hz^2));
hy0 = hy/(sqrt(hx^2+hy^2+hz^2));
hz0 = hz/(sqrt(hx^2+hy^2+hz^2));

%arranging arrays for plotting
Hx(:,k) = hx;
Hy(:,k) = hy;
Hz(:,k) = hz;
Hx0(:,k) = hx0;
Hy0(:,k) = hy0;
Hz0(:,k) = hz0;
end

t = 0:dt:(N-1)*dt;%constructing time axis

figure(1);
plot(t,Hx);
grid on;
title('Variation of x-component of the Magnetic Field Vector of Earth with respect to time');
xlabel('Time (s)');
ylabel('H_{x} (Wb/m^{2})');

figure(2);
plot(t,Hy);
grid on;
title('Variation y-component of the Magnetic Field Vector of Earth with respect to time');
xlabel('Time (s)');
ylabel('H_{y} (Wb/m^{2})');

figure(3);
plot(t,Hz);
grid on;
title('Variation z-component of the Magnetic Field Vector of Earth with respect to time');
xlabel('Time (s)');
ylabel('H_{z} (Wb/m^{2})');

figure(4);
plot(t,Hx0);
grid on;
title('Variation x-component of the Direction Cosine with respect to time');
xlabel('Time (s)');
ylabel('H_{x_{0}}');

figure(5);
plot(t,Hy0);
grid on;
title('Variation y-component of the Direction Cosine with respect to time');
xlabel('Time (s)');
ylabel('H_{y_{0}}');

figure(6);
plot(t,Hz0);
grid on;
title('Variation z-component of the Direction Cosine with respect to time');
xlabel('Time (s)');
ylabel('H_{z_{0}}');
