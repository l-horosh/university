%% Newton Method 2D
clear
clc

clc;
f = @(x) x.^3-6*x.^2+4.*x+11;
grad2 = @(x) 3*x.^2-12*x+4;
hess2 = @(x) 6*x-12;
x0 = 5;
z = newton(f, grad2, hess2, x0);

t = linspace(0, 6, 100);
plot(t, f(t),'m', z, f(z), 'ro');
xlabel('x')
ylabel('f(x)')