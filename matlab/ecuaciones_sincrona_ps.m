function dx=ecuaciones_sincrona_ps(t,X)
global voltajes corrientes tiempo enlaces velocidad par_elec angulo  V R Vm Linv tm wb h p q tipo_falla t_falla duracion_falla tm_falla;
dx=zeros(8,1);

Lambda=X(1:6)';  %% '


G=[X(2);-X(1);0;0;0;0];
w=X(7);
delta=X(8);
I=Linv*Lambda;
te=(Lambda(1)*I(2)-Lambda(2)*I(1));

if tipo_falla~=0&&t>=t_falla&&t<=(t_falla+duracion_falla/60)
    switch tipo_falla
        case 1
            V(1)=0;
            V(2)=0;
            tmn=tm;
        case 2
            tmn=tm_falla;
            V(1)=Vm*sin(delta);
            V(2)=Vm*cos(delta);
    end
else
    V(1)=Vm*sin(delta);
    V(2)=Vm*cos(delta);
    tmn=tm;
end

if q==1 %nuevas variables deseadas
    %actualizción de nuevos voltajes
    voltajes(p,:)=V';
    %nuevos enlaces
    enlaces(p,:)=Lambda';
    velocidad(p)=w;
    angulo(p)=delta;
    corrientes(p,:)=I';
    par_elec(p)=te;    %ecuación algebraica del par eléctrico
    tiempo(p)=t;    
end

dx(1:6)=(V-R*I+(w/wb)*G)*wb;         %ecuaciones diferenciales de los enlaces de flujo
dx(7)=(1/(2*h))*(tmn-te)*wb;       %ecuacion de oscilacion
dx(8)=w-wb;                       %calculo del angulo

if q==4;
    p=p+1;
    q=1;
else
    q=q+1;
end

