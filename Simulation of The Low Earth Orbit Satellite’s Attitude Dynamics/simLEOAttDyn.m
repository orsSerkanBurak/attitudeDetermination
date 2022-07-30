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
dt = 0.1; %the sample time(s)
N = 54000; %iteration number
W = zeros(3,N); %pre-allocation
for k = 1:N
%the angular velocities
w(1) = w(1) + ((dt/Jx)*(Jy - Jz)*w(3)*w(2)) + ((dt/Jx)*N_T);
w(2) = w(2) + ((dt/Jy)*(Jz - Jx)*w(1)*w(3)) + ((dt/Jy)*N_T);
w(3) = w(3) + ((dt/Jz)*(Jx - Jy)*w(1)*w(2)) + ((dt/Jz)*N_T);
    
%the euler angles
phi = phi + dt*((((w(2)*sin(phi)) + (w(3)*cos(phi)))*tan(theta)) + w(1));
theta = theta + dt*((w(2)*cos(phi)) - (w(3)*sin(phi)) + (w_orbit));
psi = psi + dt*(((w(2)*sin(phi)) + (w(3)*cos(phi)))*sec(theta));

%transformation matrix A
A = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);((-cos(phi)*sin(psi)) + (sin(phi)*sin(theta)*cos(psi))) ((cos(phi)*cos(psi)) + (sin(phi)*sin(theta)*sin(psi))) (sin(phi)*cos(theta));((sin(phi)*sin(psi))+(cos(phi)*sin(theta)*cos(psi))) ((-sin(phi)*cos(psi))+(cos(phi)*sin(theta)*sin(psi))) (cos(phi)*cos(theta))];

%arranging arrays for plotting
W(:,k) = w;
Phi(k) = phi;
Theta(k) = theta;
Psi(k) = psi;
end

t = 0:dt:(N-1)*dt;

figure(1);
plot(t,W(1,:));
grid on;
title('Angular Velocity x-axis');
xlabel('Time (s)');
ylabel('\omega_{x} (rad/s)');

figure(2);
plot(t,W(2,:));
grid on;
title('Angular Velocity y-axis');
xlabel('Time (s)');
ylabel('\omega_{y} (rad/s)');

figure(3);
plot(t,W(3,:));
grid on;
title('Angular Velocity z-axis');
xlabel('Time (s)');
ylabel('\omega_{z} (rad/s)');

figure(4)
plot(t,Phi); 
title('Roll Angle'); 
xlabel('Time (second)');
ylabel('Roll angle, \phi (rad)');
grid on;

figure(5)
plot(t,Theta); 
title('Pitch Angle'); 
xlabel('Time (second)');
ylabel('Pitch angle, \theta (rad)');
grid on;

figure(6)
plot(t,Psi); 
title('Yaw Angle'); 
xlabel('Time (second)');
ylabel('Yaw angle, \psi (rad)');
grid on;