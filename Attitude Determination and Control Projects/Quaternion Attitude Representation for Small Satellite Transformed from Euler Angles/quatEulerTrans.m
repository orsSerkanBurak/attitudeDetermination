clc;clear, close all;

n = 20; %sequence number

%initial data of the attitude angles(rad)
phi = -0.01 - 0.002 * n;
theta = 0.01 + 0.002 * n;
psi = 0.005 + 0.002 * n;

%initial data of the angular velocities of the satellite
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
q = [q1;q2;q3;q4];

%quaternion rates
q1Dot = -0.5*((q2*w(1)) + (q3*w(2)) + (q4*w(3)));
q2Dot = 0.5*((q1*w(1)) - (q4*w(2)) + (q3*w(3)));
q3Dot = 0.5*((q4*w(1)) + (q1*w(2)) - (q2*w(3)));
q4Dot = -0.5*((q3*w(1)) - (q2*w(2)) - (q1*w(3)));
qDot = [q1Dot;q2Dot;q3Dot;q4Dot];

Q(:,k) = q;
Qdot (:,k) = qDot;
end

t = 0:dt:(N-1)*dt; %constructing time axis

%plotting before applying nearest neighboring method
figure(1);
plot(t,Q(1,:)); 
title('Quaternion 1'); 
xlabel('Time(second)');
ylabel('q_{1}');
grid on;

figure(2);
plot(t,Q(2,:)); 
title('Quaternion 2'); 
xlabel('Time(second)');
ylabel('q_{2}');
grid on;

figure(3);
plot(t,Q(3,:)); 
title('Quaternion 3'); 
xlabel('Time(second)');
ylabel('q_{3}');
grid on;

figure(4);
plot(t,Q(4,:)); 
title('Quaternion 4'); 
xlabel('Time(second)');
ylabel('q_{4}');
grid on;

figure(5);
plot(t,Qdot(1,:)); 
title('Rate of Quaternion 1'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{1}$','interpreter','latex');
grid on;

figure(6);
plot(t,Qdot(2,:)); 
title('Rate of Quaternion 2'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{2}$','interpreter','latex');
grid on;

figure(7);
plot(t,Qdot(3,:)); 
title('Rate of Quaternion 3'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{3}$','interpreter','latex');
grid on;

figure(8);
plot(t,Qdot(4,:)); 
title('Rate of Quaternion 4'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{4}$','interpreter','latex');
grid on;

%nearest neighboring method
for j = 1:1:4
    for i = 1:N-1
        if (Q(j,i)*Q(j,i+1))<0
            Q(j,i+1) = -1*Q(j,i+1);
        else
            Q(j,i+1) = Q(j,i+1);
        end
    end
end

%plotting after applying nearest neighboring method

figure(9);
plot(t,Q(1,:)); 
title('Quaternion 1'); 
xlabel('Time(second)');
ylabel('q_{1}');
grid on;

figure(10);
plot(t,Q(2,:)); 
title('Quaternion 2'); 
xlabel('Time(second)');
ylabel('q_{2}');
grid on;

figure(11);
plot(t,Q(3,:)); 
title('Quaternion 3'); 
xlabel('Time(second)');
ylabel('q_{3}');
grid on;

figure(12);
plot(t,Q(4,:)); 
title('Quaternion 4'); 
xlabel('Time(second)');
ylabel('q_{4}');
grid on;

figure(13);
plot(t,Qdot(1,:)); 
title('Rate of Quaternion 1'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{1}$','interpreter','latex');
grid on;

figure(14);
plot(t,Qdot(2,:)); 
title('Rate of Quaternion 2'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{2}$','interpreter','latex');
grid on;

figure(15);
plot(t,Qdot(3,:)); 
title('Rate of Quaternion 3'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{3}$','interpreter','latex');
grid on;

figure(16);
plot(t,Qdot(4,:)); 
title('Rate of Quaternion 4'); 
xlabel('Time(second)');
ylabel('$\dot{q}_{4}$','interpreter','latex');
grid on;