
n = 100;
rng('default')
c = (-2 + 2*randn(n, 1));
f = @(u) norm(u)^2 - 3/4;
rho = supportLebesgue(f, optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp'), n); % ищем опорную функцию для множества

a = -1 + 2*randn(n, 1);
a = a/norm(a);
b = -1 + 2*randn(n, 1);
%b = b - (a.'*b) .* a;
b = b/norm(b);

%% %% Conditional Gradients
%% пример 1
rng default
A = (-2 + 2*randn(n, n));
J = @(u) norm(A*u)^2 + norm(u-a)^2 + norm(u-b)^2 + c'*u;
dJ = @(u) 2*(A')*A*u + (u-a) + (u-b) + c;
u0 = b; % начальное приближение

[u, J] = cond_grad(rho, u0, 0.0001, J, dJ);
disp(u)
disp("J = ")
disp(J)

%% пример 2
rng default
A = (-2 + 2*randn(n, n));
J = @(u) c'*u + u' * A * u;
dJ = @(u) c + (A + A') * u;
u0 = a;

[u, J] = cond_grad(rho, u0, 0.0001, J, dJ);
disp(u)
disp("J = ")
disp(J)

%% Newton's Method
%% пример
rng default
A = (-2 + 2*randn(n, n));
f = (-2 + 2*randn(n, 1));
J = @(u) 0.5*u'*A*u - 2*f'*u;
dJ = @(u) A'*A*u-f;
d2J = @(u) A'*A;

u0 = 3*f;
u = newton(dJ, d2J, u0);
disp(u)
disp("J = ")
disp(J(u))

%% Simplex Method

A = [1 0 4 1 -1; 1 1 6 2 0; 0 0 1 10 0];
b = [1; -1; 1];
c = [1; 1; -1; 0; 1];
z0 = [3; 0; 11; -1; 9];
[u, J] = simplex(A, b, c, z0);
disp(u)
disp("J = ")
disp(J)

%%

A = [1 0 0 1 -1; 1 1 0 2 0; 0 0 1 1 0];
b = [1; 3; 1];
v = [3 0 1 0 2]';
c = [1 1 -1 0 1]';
[u, J] = simplex(A, b, c, v);
disp(u)
disp("J = ")
disp(J)
