v1B = [0.8273, 0.5541, -0.0920];
v2B = [-0.8285, 0.5522, -0.0955];
v1N = [-0.1517, -0.9669, 0.2050];
v2N = [-0.8393, 0.4494, -0.3044];

t1B = v1B;
t2B = cross(v1B, v2B);
t2B = t2B / norm(t2B);
t3B = cross(t1B, t2B);

t1N = v1N;
t2N = cross(v1N, v2N);
t2N = t2N / norm(t2N);
t3N = cross(t1N, t2N);

BbarT = [t1B; t2B; t3B];
NT = [t1N; t2N; t3N];
BbarN = BbarT * NT';

% Displaying the results
disp('t1B:');
disp(t1B);
disp('t2B:');
disp(t2B);
disp('t3B:');
disp(t3B);
disp('t1N:');
disp(t1N);
disp('t2N:');
disp(t2N);
disp('t3N:');
disp(t3N);
disp('BbarT:');
disp(BbarT);
disp('BbarN:');
disp(BbarN);
