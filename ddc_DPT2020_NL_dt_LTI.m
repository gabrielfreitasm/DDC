clear
close all
clc

warning off

% Illustrative example: SSF NL system

% inverted pendulum discretized by Euler method
n = 2; % two states
m = 1; % single input

Tmin = (m + 1)*n + m;
T = 5; % simulation time

Ts = 0.1; % sampling time

mp = 1;
l = 1;
g = 9.8;
muu = 0.01;

u = -0.1 + 0.2*rand(m, T); % random input sequence of length T within [-0.1, 0.1]

u0 = u(:, 1);
x10 = -0.1 + 0.2*rand; % random initial conditions state x1 within [-0.1, 0.1]
x20 = -0.1 + 0.2*rand; % random initial conditions state x2 within [-0.1, 0.1]

x1(1) = x10 + Ts*x20;
x2(1) = Ts*g/l*sin(x10) + (1 - Ts*muu/(mp*l^2))*x20 + Ts/(mp*l^2)*u0;
for i = 1:T-1
    x1(i+1) = x1(i) + Ts*x2(i);
    x2(i+1) = Ts*g/l*sin(x1(i)) + (1 - Ts*muu/(mp*l^2))*x2(i) + Ts/(mp*l^2)*u(i);
end

d(:, 1) = [0; Ts*g/l*(sin(x10) - x10)];
for i = 2:T
    d(:, i) = [0; Ts*g/l*(sin(x1(i)) - x1(i))];
end

X0 = [x10, x1(1:T-1); x20, x2(1:T-1)]; % state samples
U0 = [u0, u(1:T-1)]; % input samples
X1 = [x1; x2]; % dynamics samples
D0 = d; % higher-order terms

% assumption 5:
cvx_begin sdp
    
    variable gam

    minimize gam

    subject to
    gam >= 0;
    (D0*D0') <= gam*(X1*X1') ;

cvx_end

% eq. (53)
% alph = 0.0422;
cvx_begin sdp

    variable Q(T, n)
    variable alph

    maximize alph
    
    subject to
    alph >= 0;

    [X0*Q - alph*(X1*X1'), X1*Q; Q'*X1', X0*Q] >= 0;
    
    [eye(T), Q; Q', X0*Q] >= 0;

cvx_end

% test condition gam < alph^2/(4 + 2*alph)
disp('Max SNR condition satisfied? (conservative condition)')
if gam < alph^2/(4 + 2*alph)
    disp('Yes')
else
    disp('No')
end

% eq. (53)
Gk = Q/(X0*Q);
K = U0*Gk;

% test closed-loop stability
A = [1, 0.1; 0.98, 0.999];
B = [0; 0.1];

abs(eig(A + B*K))

return