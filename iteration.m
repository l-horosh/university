% Function for Newton Method
function u_opt = iteration(A, f, u0)
    epsilon=0.001;
    n = length(f);
    J = @(u) 0.5*u'*A*u - f'*u; %функционал
    dJ = @(u) A*u-f; %Производная функционала
    scalar = @(a,b) a'*b; %скалярное произведение
    norma = @(a) sqrt(scalar(a,a)); %норма вектора
    alpha_formula = @(u,p) -scalar(dJ(u),p)/scalar(A*p,p); %формула (2) со страницы 38 лекций 
    beta_formula = @(u,p) scalar(dJ(u),A*p)/scalar(A*p,p); %формула (4\beta) со страницы 39 лекций
    
    u_tmp = zeros(n,2); %в первом столбце хрантся u_k, а во втором u_{k+1}
    p_tmp = zeros(n,2); %в первом столбце храним p_k, а во втором p_{k+1}
    u_tmp(:,1) = u0;
    p_tmp(:,1) = -dJ(u0);

    while norma(dJ(u_tmp(:,1)))/norma(f) > epsilon %проверяем размер невязки
        alpha = alpha_formula(u_tmp(:,1),p_tmp(:,1)); %вычисялем альфа
        u_tmp(:,2) = u_tmp(:,1)+alpha*p_tmp(:,1); %вычисляем очередное uk по формуле (4u) со страницы 39 лекций
        beta = beta_formula(u_tmp(:,2),p_tmp(:,1)); %вычисляем бета
        p_tmp(:,2) = -dJ(u_tmp(:,2)) + beta*p_tmp(:,1); %вычисляем очередное uk по формуле (4p) со страницы 39 лекций
        u_tmp = circshift(u_tmp,[0,-1]); %циклическое смещение матрицы влево на 1
        p_tmp = circshift(p_tmp,[0,-1]); %циклическое смещение матрицы влево на 1
    end
    u_opt=u_tmp(:,1);
end