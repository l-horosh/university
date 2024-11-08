x0=[0.0092;-0.9];
sigma=0.3;
dzeta=0.64;
alpha=-0.5; beta=0.5; gamma=1; r=1;
A=[0.5,2.0;-0.5,0]; B=[1, 0;0, 1]; f=[0;0];
dfi=0.01;   %  шаг по углу начального условия для psi
NX1=1000;   %  число разбиений контура целевой области
Tmax=3;     %  максимальное время моделирования траектории
Rtol=1e-9;  %  относительная точность при вычислении траектории






