% Function for Conjugate gradients Method
% J(u) = ||Au-f|| + <c,u> -> inf
function [u_min, J_min] = conj_grad(A, f, u0, eps)
J = @(x) sum((A*x - f).^2);
grad = @(x) 2*A'*(A*x - f);
n_max = 1000;
dim = numel(u0);
u = u0;
n = 0;

p = -grad(u0);
while (sum(grad(u).^2) >= eps^2) && (n < n_max)
     u_prev = u;
     alpha = - (grad(u_prev)' * p) / (2 * p' * A * p);
     u = u_prev + alpha * p;
     gu = grad(u);
     beta = (gu' * A * p) / (p' * A * p);
     p = -gu + beta * p;
     n = n + 1;
     if dim == 2
         subplot(1,1,1);
         plot([u_prev(1) u(1)], [u_prev(2) u(2)],'o:r');
         hold on;
     end
end

if n == n_max
    disp('Algorithm got into a cycle');
    u_min = [];
    J_min = [];
    return;
else
    u_min = u;
    J_min = J(u);
    disp(['Number of iterations: ', num2str(n)]);
    disp('u* = ');
    disp(u_min);
    disp('J* = ');
    disp(J_min);
    return;
end
end