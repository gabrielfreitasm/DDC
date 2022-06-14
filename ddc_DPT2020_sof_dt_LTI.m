clear
close all
clc

warning off

% Illustrative example: SOF

Ts = 1; % sampling time

a1 = 1;
a2 = - 2.311;
a3 = 2.623;
a4 = - 2.311;

a = [a1; a2; a3; a4];

aux = [zeros(1, size(a,1)-1);
    eye(size(a,1)-1)];

b1 = 0.039;
b2 = 0.383;
b3 = 0.383;
b4 = 0.039;

b = [b1; b2; b3; b4];

A = [[aux'; -a'], [zeros(size(a,1)-1, size(b,1)); b'];
    zeros(size(a,1), size(a,1)), [aux'; zeros(1, size(b,1))]];

B = [zeros(2*size(a,1)-1,1);
    1];

C = [-a', b'];

n = size(A,1);
m = size(B,2);

Tmin = 2*(n/2) + 1;

T = Tmin;

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

cvx_begin sdp

%     variable P(size(A,1), size(A,2)) semidefinite
    variable Q(T, size(A,2))

    % eq. (74)
    [X0*Q, X1*Q; Q'*X1', X0*Q] >= 1e-5*eye(2*size(Q,2));

    % numerical implementation
%     [P, X1*Q; Q'*X1', P] >= 1e-5*eye(2*size(P));
%     X0*Q == P;

cvx_end

% eq. (74)
Gk = Q/(X0*Q);
K = U0*Gk;

% numerical implementation
% Gk = Q/P;
% K = U0*Gk;

abs(eig(A + B*K))

return