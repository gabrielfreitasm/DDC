clear
close all
clc

warning off

% Illustrative example: SSF with noise

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

w = -0.01*ones(n, T+1) + 0.02*rand(n, T+1); % noise random sequence of amplitude 0.02

w0 = w(:, 1);
u0 = u(:, 1);
x0 = rand(n, 1); % random initial conditions

x(:, 1) = A*(x0 + w0) + B*u0;

for i = 1:T-1
    if i == 1
        x(:, i+1) = A*x(:, i) + B*u(:,i);
    else
        x(:, i+1) = A*(x(:, i) + w(:, i)) + B*u(:,i);
    end
end

X0 = [x0, x(:, 1:T-1)]; % state samples
W0 = w(:, 1:T);
U0 = [u0, u(:, 1:T-1)]; % input samples
X1 = [B, A]*[U0; X0]; % dynamics samples
W1 = w(:, 2:T+1);
Z0 = X0 + W0;
Z1 = X1 + W1;

% eq. (35)
cvx_begin sdp

    variable Q(T, size(A,2))
    variable alph
    
    maximize alph

    subject to
    
    alph >= 0;

    [Z0*Q - alph*(Z1*Z1'), Z1*Q; Q'*Z1', Z0*Q] >= 0;
    
    [eye(T), Q; Q', Z0*Q] >= 0;

cvx_end

% assumption 2: max SNR gamma
R0 = A*W0 - W1;
cvx_begin sdp quiet
    
    variable gam
    
    minimize gam
    
    subject to
    gam >= 0;
    (R0*R0') <= gam*(Z1*Z1') ;

cvx_end

% test condition gam < alph^2/(4 + 2*alph)
disp('Max SNR condition satisfied? (conservative condition)')
if gam < alph^2/(4 + 2*alph)
    disp('Yes')
else
    disp('No')
end

% eq. (35)
Gk = Q/(Z0*Q);
K = U0*Gk;

abs(eig(A + B*K))

return