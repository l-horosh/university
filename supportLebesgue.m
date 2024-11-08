function res = supportLebesgue(f, opts, n)
    res = @rho;
    function [val, point] = rho(y)
        nonlcon = @(x) const(f(x));
        maxim = fmincon(@(x) -(y' * x), zeros(n, 1), [], [], [], [], [], [], nonlcon, opts); %т.к. ищем супремум, то в fmin находим -inf
        val = y' * maxim;
        point = maxim;
    end
end