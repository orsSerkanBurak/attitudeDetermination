clc;clear;close all;

n = 20; %sequence number

%initial data of the attitude angles(rad)
phi = -0.01 - 0.002 * n;
theta = 0.01 + 0.002 * n;
psi = 0.005 + 0.002 * n;

%initial data of the angular velocities of satellite
wX = 0.0002 + 0.00001 * n;
wY = 0.0003 + 0.00001 * n;
wZ = 0.0004 + 0.00001 * n;
w = [wX;wY;wZ];

%initial moments of inertia of the satellite(m^4)
Jx = 2.1*(10^(-3));
Jy = 2*(10^(-3));
Jz = 1.9*(10^(-3));

w_orbit = 0.0011; %The angular orbit velocity of satellite(rad/s)
N_T = 3.6*(10^(-10)); %The disturbance torque acting on the satellite(Nm)

Me = 7.943e15; %the magnetic dipole moment of the Earth in Wb.m
mu = 3.98601e14; %the Earth Gravitational constant in m^3/s^2
i = (80 + (0.5*n))*(pi/180); %inclination converted from deg to rad
we = 7.29e-5; %the spin rate of the Earth in rad/s
epsilon = 11.7*(pi/180); %the magnetic dipole tilt converted to rad
r0 = (6378.14 + 500 + (2*n))*1000; %the distance between the center of mass of the satellite and the Earth in m
w0 = sqrt(mu/r0^3); %the angular velocity of the orbit with respect to the inertial frame in rad/s

dt = 0.1; %the sample time(s)
N = 54000; %iteration number

% pre-allocation
Bm0=zeros(3,N);
Bm0_star=zeros(3,N);
Bm=zeros(3,N);
Bm_star=zeros(3,N);

for k = 1:N
% the angular velocities
w(1) = w(1) + ((dt/Jx)*(Jy - Jz)*w(3)*w(2)) + ((dt/Jx)*N_T);
w(2) = w(2) + ((dt/Jy)*(Jz - Jx)*w(1)*w(3)) + ((dt/Jy)*N_T);
w(3) = w(3) + ((dt/Jz)*(Jx - Jy)*w(1)*w(2)) + ((dt/Jz)*N_T);
    
% the euler angles
phi = phi + dt*((((w(2)*sin(phi)) + (w(3)*cos(phi)))*tan(theta)) + w(1));
theta = theta + dt*((w(2)*cos(phi)) - (w(3)*sin(phi)) + (w_orbit));
psi = psi + dt*(((w(2)*sin(phi)) + (w(3)*cos(phi)))*sec(theta));

% transformation matrix A
A = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);((-cos(phi)*sin(psi)) + (sin(phi)*sin(theta)*cos(psi))) ((cos(phi)*cos(psi)) + (sin(phi)*sin(theta)*sin(psi))) (sin(phi)*cos(theta));((sin(phi)*sin(psi))+(cos(phi)*sin(theta)*cos(psi))) ((-sin(phi)*cos(psi))+(cos(phi)*sin(theta)*sin(psi))) (cos(phi)*cos(theta))];

% components of Earth Magnetic Field
hx = (Me/r0^3)*((cos(w0*k)*((cos(epsilon)*sin(i))-(sin(epsilon)*cos(i)*cos(we*k)))-(sin(w0*k)*sin(epsilon)*sin(we*k))));
hy = (-(Me/r0^3))*((cos(epsilon)*cos(i))+(sin(epsilon)*sin(i)*cos(we*k)));
hz = ((2*Me)/r0^3)*((sin(w0*k)*((cos(epsilon)*sin(i))-(sin(epsilon)*cos(i)*cos(we*k)))-(2*(sin(w0*k)*sin(epsilon)*sin(we*k)))));
h =[hx;hy;hz]*10^6; %converted to microTesla

% components of Direction Cosine Matrix
hx0 = hx/(sqrt(hx^2+hy^2+hz^2));
hy0 = hy/(sqrt(hx^2+hy^2+hz^2));
hz0 = hz/(sqrt(hx^2+hy^2+hz^2));
h0=[hx0;hy0;hz0];

% standard deviation of each magnetometer errors
sigma_mC = 0.008; % in terms of direction cosine
sigma_m = 1.66; % in terms of projectiles

% magnetometer bias vectors
bc = [0.04;0.06;0.08]; % in terms of direction cosine
b = [3;5;6]; % in terms of projectiles

% zero mean Gaussian white noises
rng(0);
vc = [sigma_mC*randn;sigma_mC*randn;sigma_mC*randn]; %in terms of direction cosine
v = [sigma_m*randn;sigma_m*randn;sigma_m*randn];

%magnetometer measurements in terms of projectiles
Bm(:,k) = ((A*h)+b+v); % with bias and noise
Bm_star(:,k) = (A*h); % without bias and noise

%magnetometer measurements in terms of direction cosines
Bm0(:,k)=(A*h0)+bc+vc; % with bias and noise
Bm0_star(:,k)=(A*h0); % without bias and noise
end

t = 0:dt:(N-1)*dt; %constructing time axis

figure(1)
bmx = plot(t,Bm(1,:));
grid on
hold on
bmsx = plot(t,Bm_star(1,:));
hold off
title('Variation x-component of the Magnetometer Measurement in terms of Projectile with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{x}} (microTesla, \muT)');
legend([bmx,bmsx],'with noise, B_{m}','without noise,B_{m}^{*}');

figure(2)
bmy = plot(t,Bm(2,:));
grid on
hold on
bmsy = plot(t,Bm_star(2,:));
hold off
title('Variation y-component of the Magnetometer Measurement in terms of Projectiles with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{y}}(microTesla, \muT)');
legend([bmy,bmsy],'with noise, B_{m}','without noise,B_{m}^{*}');

figure(3)
bmz = plot(t,Bm(3,:));
grid on
hold on
bmsz = plot(t,Bm_star(3,:));
hold off
title('Variation z-component of the Magnetometer Measurement in terms of Projectiles with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{z}}(microTesla, \muT)');
legend([bmz,bmsz],'with noise, B_{m}','without noise,B_{m}^{*}');

figure(4)
bm0x = plot(t,Bm0(1,:));
grid on
hold on
bm0sx = plot(t,Bm0_star(1,:));
hold off
title('Variation x-component of the Magnetometer Measurement in terms of Direction Cosines with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{0_{x}}}');
legend([bm0x,bm0sx],'with noise, B_{m_{0}}','without noise,B_{m_{0}}^{*}');

figure(5)
bm0y = plot(t,Bm0(2,:));
grid on
hold on
bm0sy = plot(t,Bm0_star(2,:));
hold off
title('Variation y-component of the Magnetometer Measurement in terms of Direction Cosines with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{0_{y}}}');
legend([bm0y,bm0sy],'with noise, B_{m_{0}}','without noise,B_{m_{0}}^{*}');

figure(6)
bm0z = plot(t,Bm0(3,:));
grid on
hold on
bm0sz = plot(t,Bm0_star(3,:));
hold off
title('Variation z-component of the Magnetometer Measurement in terms of Direction Cosines with respect to time');
xlabel('Time (s)');
ylabel('B_{m_{0_{z}}}');
legend([bm0z,bm0sz],'with noise, B_{m_{0}}','without noise,B_{m_{0}}^{*}');
