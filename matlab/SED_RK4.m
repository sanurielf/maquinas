function [t,y] = SED_RK4(f,rango,y0,h)
%Metodo de runge kutta para resolver un sistema de ecuaciones diferenciales expresando en la forma vectorial y’(t) = f(t,y(t))
% para un rango = [t0,tf] con condiciones iniciales yo=[y1(i) y2(i)...]. N=
% numero de muestras.


% proteccion contra N negativa o en caso de que no se escoja N se asgina un
% tamapo de paso de .001
if nargin < 4 || h <= 0
    h=.00001;
end
N=(rango(2)-rango(1))/h;
if nargin < 3
    y0 = 0;
end

y(1,:) = y0(:)'; %hace un vector con las condiciones iniciales.

h = (rango(2) - rango(1))/N;
t =rango(1)+[0:N]'*h;% vector tiempo

for k = 1:N
    f1 = h*feval(f,t(k),y(k,:));                f1=f1(:)';
    f2 = h*feval(f,t(k) + h/2,y(k,:) + f1/2);   f2=f2(:)';
    f3 = h*feval(f,t(k) + h/2,y(k,:) + f2/2);   f3=f3(:)';
    f4 = h*feval(f,t(k) + h,y(k,:) + f3);       f4=f4(:)';
    y(k + 1,:) = y(k,:) + (f1 + 2*(f2 + f3) + f4)/6;
end