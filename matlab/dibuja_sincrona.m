function dibuja_sincrona(monitor)

load('resultados_sincrona1.mat');
switch monitor
    case 1 %monitor casa
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(2) scrsz(3) scrsz(4)])
        subplot(2,4,1)
        %voltajes
        plot(tiempo,voltajes(:,1));
        hold all;
        plot(tiempo,voltajes(:,2));
        legend('Vd','Vq');
        title('Voltajes');
        xlabel('tiempo');
        ylabel('Voltaje pu');
        %corrientes
        subplot(2,4,2)
        plot(tiempo,corrientes(:,1));
        hold all;
        plot(tiempo,corrientes(:,2));
        legend('Id','Iq');
        title('Corrientes');
        xlabel('tiempo');
        ylabel('Corriente pu');
        %velocidad
        subplot(2,4,3)
        plot(tiempo,velocidad);
        legend('Velocidad');
        title('Velocidad del rotor');
        xlabel('tiempo');
        ylabel('RPM');
        %angulo
        subplot(2,4,4)
        plot(tiempo,angulo);
        legend('Angulo delta');
        title('Angulo de carga');
        xlabel('tiempo');
        ylabel('Grados');
        %par
        subplot(2,4,5)
        plot(tiempo,par_elec);
        legend('Par eléctrico');
        title('Par eléctrico');
        xlabel('tiempo');
        ylabel('N/m pu');
        %velocidad-angulo
        subplot(2,4,6)
        plot(velocidad,angulo);
        legend('Velocidad vs Angulo');
        title('Velocidad vs Angulo');
        ylabel('Grados');
        xlabel('RPM');
        %angulo-par
        subplot(2,4,7);
        plot(angulo,par_elec);
        legend('Angulo vs Par eléctrico');
        title('Angulo vs Par eléctrico');
        xlabel('Grados');
        ylabel('N/m pu');
    case 2 %monitor laptop
        %voltajes
        plot(tiempo,voltajes(:,1));
        hold all;
        plot(tiempo,voltajes(:,2),'r');
        legend('Vd','Vq');
        title('Voltajes');
        xlabel('tiempo');
        ylabel('Voltaje pu');
        %corrientes
        figure;
        plot(tiempo,corrientes(:,1));
        hold all;
        plot(tiempo,corrientes(:,2),'r');
        legend('Id','Iq');
        title('Corrientes');
        xlabel('tiempo');
        ylabel('Corriente pu');
        %velocidad
        figure;
        plot(tiempo,velocidad);
        legend('Velocidad');
        title('Velocidad del rotor');
        xlabel('tiempo');
        ylabel('RPM');
        %angulo
        figure;
        plot(tiempo,angulo);
        legend('Angulo delta');
        title('Angulo de carga');
        xlabel('tiempo');
        ylabel('Grados');
        %par
        figure;
        plot(tiempo,par_elec);
        legend('Par eléctrico');
        title('Par eléctrico');
        xlabel('tiempo');
        ylabel('N/m pu');
        %velocidad-angulo
        figure;
        plot(velocidad,angulo);
        legend('Velocidad vs Angulo');
        title('Velocidad vs Angulo');
        ylabel('Grados');
        xlabel('RPM');
        %angulo-par
        figure;
        plot(angulo,par_elec);
        legend('Angulo vs Par eléctrico');
        title('Angulo vs Par eléctrico');
        xlabel('Grados');
        ylabel('N/m pu');
end
