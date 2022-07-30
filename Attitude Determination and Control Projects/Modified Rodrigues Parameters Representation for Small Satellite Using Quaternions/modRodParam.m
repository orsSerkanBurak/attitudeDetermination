clc;clear;close all;

n = 20; %sequence number

%initial data of the attitude angles(rad)
phi = -0.01 - 0.002 * n;
theta = 0.01 + 0.002 * n;
psi = 0.005 + 0.002 * n;

%initial data of angular velocities of the satellite
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
dt = 0.1; %the sample time(s)
N = 54000; %iteration number
for k = 1:N
%the angular velocities
w(1) = w(1) + ((dt/Jx)*(Jy - Jz)*w(3)*w(2)) + ((dt/Jx)*N_T);
w(2) = w(2) + ((dt/Jy)*(Jz - Jx)*w(1)*w(3)) + ((dt/Jy)*N_T);
w(3) = w(3) + ((dt/Jz)*(Jx - Jy)*w(1)*w(2)) + ((dt/Jz)*N_T);
    
%the euler angles
phi = phi + dt*((((w(2)*sin(phi)) + (w(3)*cos(phi)))*tan(theta)) + w(1));
theta = theta + dt*((w(2)*cos(phi)) - (w(3)*sin(phi)) + (w_orbit));
psi = psi + dt*(((w(2)*sin(phi)) + (w(3)*cos(phi)))*sec(theta));

%quaternions
q1 = cos(phi/2)*cos(theta/2)*cos(psi/2)+sin(phi/2)*sin(theta/2)*sin(psi/2);
q2 = sin(phi/2)*cos(theta/2)*cos(psi/2)-cos(phi/2)*sin(theta/2)*sin(psi/2);
q3 = cos(phi/2)*sin(theta/2)*cos(psi/2)+sin(phi/2)*cos(theta/2)*sin(psi/2);
q4 = cos(phi/2)*cos(theta/2)*sin(psi/2)+sin(phi/2)*sin(theta/2)*cos(psi/2);
q = [abs(q1);abs(q2);abs(q3);abs(q4)];

%modified Rodrigues parameters
p1 = q(1)/(1+q(4));
p2 = q(2)/(1+q(4));
p3 = q(3)/(1+q(4));
p = [p1;p2;p3];

%arranging arrays for plotting
P(:,k) = p;
end

t = 0:dt:(N-1)*dt;

figure(1);
plot(t,P(1,:)); 
title('Modified Rodrigues Parameter 1'); 
xlabel('Time(second)');
ylabel('p_{1}');
grid on;

figure(2);
plot(t,P(2,:)); 
title('Modified Rodrigues Parameter 2'); 
xlabel('Time(second)');
ylabel('p_{2}');
grid on;

figure(3);
plot(t,P(3,:)); 
title('Modified Rodrigues Parameter 3'); 
xlabel('Time(second)');
ylabel('p_{3}');
grid on;