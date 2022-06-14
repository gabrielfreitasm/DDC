clear
close all
clc

warning off

% Illustrative example: LQR stabilization

Ts = 0.1; % sampling time

AB = [1.178, 0.001, 0.511, - 0.403, 0.004, - 0.087;
    - 0.051, 0.661, - 0.011, 0.061, 0.467, 0.001;
    0.076, 0.335, 0.560, 0.382, 0.213, - 0.235;
    0, 0.335, 0.089, 0.849, 0.213, - 0.016];

A = AB(:, 1:4);
B = AB(:, 5:end);

n = size(A,1);
m = size(B,2);

Tmin = (m + 1)*n + m;

T = Tmin+1;

u = rand(m, T); % random input sequence of length T

u0 = u(:, 1);
x0 = rand(n, 1); % random initial conditions
x(:, 1) = A*x0 + B*u0;
for i = 1:T-1
    x(:, i+1) = A*x(:, i) + B*u(:,i);
end

X0 = [x0, x(:, 1:T-1)]; % state samples
U0 = [u0, u(:, 1:T-1)]; % input samples
X1 = [B, A]*[U0; X0]; % dynamics samples

Qx = eye(n);
R = eye(m);

cvx_begin sdp

    variable Q(T, n)
    variable X(m, m)

    minimize trace(Qx*X0*Q) + trace(X)
    
    subject to
    % eq. (27)
    [X, R^(0.5)*U0*Q; Q'*U0'*R^(0.5), X0*Q] >= 1e-5*eye(m + n);

    [X0*Q - eye(n), X1*Q; Q'*X1', X0*Q] >= 1e-5*eye(2*n);

cvx_end

% eq. (27)
Gk = Q/(X0*Q);
K = U0*Gk;

abs(eig(A + B*K))

return