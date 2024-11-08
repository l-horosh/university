clear; close all;
global alpha beta gamma r
par;
if exist("x0","var")==0 || any(size(x0)~=[2,1])
    fprintf("Ошибка входных данных: некорректное начальное множество\n");
    return;
end
if exist("A","var")==0 || any(size(A)~=[2,2])
    fprintf("Ошибка входных данных: некорректная матрица A\n");
    return;
end

Psi=linspace(0,2*pi*0.99,100);
poly1 = polyshape([-sigma+sqrt(dzeta)*cos(Psi)],[sqrt(dzeta)*sin(Psi)]);
poly2 = polyshape([sigma+sqrt(dzeta)*cos(Psi)],[sqrt(dzeta)*sin(Psi)]);
polyout = intersect(poly1,poly2);
if polyout.NumRegions<=0
    fprintf("Ошибка входных данных: множество управления пусто\n");
    return;
end
figure(1); plot(polyout); hold on;

P=[];
for psi=Psi
    l=[cos(psi), sin(psi)];
    [rho,v]=rhoP(l,sigma,dzeta);
    P=[P;v];
end
figure(1); plot(P(:,1),P(:,2)); hold on; title("Область управления");
xlim([-2,2]);
ylim([-2,2]);

X1=linspace(alpha-sqrt(gamma),alpha+sqrt(gamma),NX1);
XT=max(0,gamma-(X1-alpha).^2);
X2=beta+sqrt(XT);
X2a=beta-sqrt(XT);
%plot(X1,X2);plot(X1,X2a);
%plot([r,0,-r,0,r],[0,r,0,-r,0])
poly1 = polyshape([X1,flip(X1(2:end-1))],[X2,flip(X2a(2:end-1))]);
poly2 = polyshape([r,0,-r,0,r],[0,r,0,-r,0]);
polyout = intersect(poly1,poly2);
if polyout.NumRegions<=0
    fprintf("Ошибка входных данных: целевое множество пусто\n");
    return;
end
figure(2); hold on;
plot(polyout); plot([x0(1)],[x0(2)],'b.'),text(x0(1),x0(2),'X_0')
if isinterior(polyout,x0(1),x0(2))
    fprintf("Tmin=0 - исходная точка находится внутри целевой области\n");
    return;
end

Tmin=200; Hit=false;
P=[];
for psi=Psi
    l=[cos(psi), sin(psi)];
    [rho,v]=rhoP(l,sigma,dzeta);
    P=[P;v];
end
for fi=0:dfi:2*pi
   Y0=[x0; cos(fi); sin(fi)];
   [T,Y,hit]=traject(Y0,A,B,f,sigma,dzeta,Tmax, Rtol);
   figure(2); plot(Y(:,1),Y(:,2));
   figure(3); plot(Y(:,3),Y(:,4)); hold on;
   if T(end)<Tmin && hit
       Hit=true;
       Tmin=T(end);
       Ym=Y;
   end
end
if ~Hit
    fprintf("При таких данных задача управления неразрешима\n");
    return;
end
fprintf("Tmin=%g\n",Tmin);
figure(2); plot(Ym(:,1),Ym(:,2),'LineWidth',3);
figure(3); plot(Ym(:,3),Ym(:,4),'LineWidth',3);axis equal;
Nm=length(Ym(:,3));
u=zeros(Nm,2);
for i=1:Nm
  [rho,u(i,:)]=rhoP((B'*Ym(i,3:4)')',sigma,dzeta);
end
figure(1); plot(u(:,1),u(:,2),'r.');axis equal;

figure(2); plot([Ym(end,1),Ym(end,1)-Ym(end,3)],[Ym(end,2),Ym(end,2)-Ym(end,4)])
P=[];
for fi=0:dfi:2*pi
    l=[cos(fi), sin(fi)];
    [rho,v]=rhoX(l,polyout);
    P=[P;v];
end
plot(P(:,1),P(:,2));
axis equal;
psi1=[Ym(end,3),Ym(end,4)];
psi1=psi1/norm(psi1);
[rhoX1,svX1]=rhoX(-psi1,polyout);
%delta=abs(-psi1*[Ym(end,1);Ym(end,2)]-rhoX(-psi1,polyout));
delta=abs(-psi1*[Ym(end,1);Ym(end,2)]-rhoX1);
fprintf("Погрешность трансверсальности в конечной точке %g\n",delta);


function [Y,T,hit]=traject(Y0,A,B,f,sigma,dzeta,Tmax, Rtol)
odefun=@(t,y) rfun(y,A,B,f,sigma,dzeta);
opts = odeset('RelTol',Rtol,'Events',@inX1,'MaxStep',0.01);
[Y,T,te,ye,ie]=ode45(odefun,[0,Tmax],Y0,opts);
hit=false;
if (~isempty(ie)), hit=true; end
end

function [position,isterminal,direction] = inX1(t,y)
global alpha beta gamma r
  position = ((y(1)-alpha)^2+(y(2)-beta)^2>gamma) + ...
             (abs(y(1))+abs(y(2))>r);
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

function dy=rfun(y,A,B,f,sigma,dzeta)
dy=zeros(4,1);
dy(3:4)=-A'*y(3:4);
[rho,u]=rhoP((B'*y(3:4))',sigma,dzeta);
dy(1:2)=A*y(1:2)+B*u'+f;
end


function [rho,v]=rhoP(l,sigma,dzeta)
l=l/norm(l);
rho1=l*[-sigma,0]'+sqrt(dzeta);
rho2=l*[sigma,0]'+sqrt(dzeta);
rho=rho1; v=[-sigma,0]+l*sqrt(dzeta);
if rho>rho2, rho=rho2; v=[sigma,0]+l*sqrt(dzeta); end
yi=sqrt(dzeta-sigma^2);
if abs(v(2))>yi, rho=l*[0;yi*sign(l(2))]; v=[0,yi*sign(l(2))]; end
end

function [rho,sv]=rhoX(l,polyX)
l=l/norm(l);
rho=-1e-10;
for v=polyX.Vertices'
    if l*v >rho
      rho=l*v;
      sv=v';
    end
end    
end
