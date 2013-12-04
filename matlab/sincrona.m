function sincrona
% máquina sincrona de polos salientes en el estator, un devanado de campo  y
% dos devanados de amortiguamiento s y t colocados sobre el eje d y q
% Autor: Uriel Sandoval
% Dic-2011

clear all;
clc;
global   ti tf pp V voltajes corrientes tiempo enlaces velocidad par_elec angulo p q variables_estado tipo_falla t_falla duracion_falla tm_falla; %variables globales

ti=1;
tf=10;
paso=.001;
maquina=2;
condiciones=1;
%Definición de parámetros
param_gensal(maquina,condiciones);%carga parametros de la máquina y condiciones iniciales

% *************************************** modulo de fallas
tipo_falla=2; % 1=falla trifásica, 2=cambio de par mecánico
t_falla=2; %segundos
duracion_falla=6; %ciclos
tm_falla=0; %aplica una falla a la máquina en un determinado tiempo con un duración en ciclos

p=1;
q=1;
%Resolver conjunto de ecuaciones diferenciales
SED_RK4('ecuaciones_sincrona_ps',[ti tf],variables_estado,paso);
velocidad=120*velocidad/(2*pi*pp); %conversion de rad/seg a rpm
angulo=angulo*180/pi;                   %conversion de rad a grados
save('resultados_sincrona1.mat','voltajes','corrientes','tiempo','enlaces','velocidad','par_elec','angulo');