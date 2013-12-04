function param_gensal(modelo,condiciones)
global ws wb pp h V tm Vm R L Linv variables_estado

switch modelo
    case 1
        frec=60;                                                                                                    %frecuencia del sistema;
        wb=2*pi*frec;
        ws=wb;      %velocidad base
        pp=2;      %pares de polos
        h=.3568;                                                                                                        %constante de inercia
        
        tm=.85;                                                                                                        %par mecánico
        Vm=1;                                                                                                       %voltaje máximo
        
        %Parametros fisicos de los devanados de la máquina (resistencias e
        %inductancias
        r=8.979e-3;
        rf=.00206;
        rs=.2826;
        rt=.02545;
        
        
        ll=.05;
        lmd=2.35;
        lmq=1.72;
        llfd=.511;
        llkd=3.738;
        llkq=.2392;
        
        %inductancias propias
        ld=ll+lmd;
        lq=ll+lmq;
        l0=ll;
        lf=llfd+lmd;
        ls=llkd+lmd;
        lt=llkq+lmq;
        %inductancias mutuas
        mdf=lmd;
        mds=lmd;
        mqt=lmq;
        switch condiciones
            case 1
                
                V=[0                                                                                       %vector de voltajes inicial
                    0
                    0
                    .001;
                    0
                    0];
                enlaces=[0 0 0 0 0 0];
                velocidad=wb;
                angulo=0;
        end
    case 2
        frec=60;                                                                                                    %frecuencia del sistema;
        wb=2*pi*frec;
        ws=wb;      %velocidad base
        pp=32;      %pares de polos
        h=.5*7.5;                                                                                                        %constante de inercia
        tm=.8519;                                                                                                        %par mecánico
        Vm=1;                                                                                                       %voltaje máximo
        %Parametros fisicos de los devanados de la máquina (resistencias e
        %inductancias
        r=.0019;
        rf=4.1e-4;
        rs=1.41e-2;
        rt=1.36e-2;
        
        ll=.12;
        llfd=.2049;
        llkd=.16;
        llkq=.1029;
        
        %inductancias propias
        ld=.85;
        lq=.48;
        l0=.12;
        lmd=ld-ll;
        lmq=lq-ll;
        lf=llfd+lmd;
        ls=llkd+lmd;
        lt=llkq+lmq;
        %inductancias mutuas
        mdf=lmd;
        mds=lmd;
        mqt=lmq;
        switch condiciones
            case 1
                V=[.3086                                                                                       %vector de voltajes inicial
                    .9512
                    0
                    8.9936e-4;
                    0
                    0];
                enlaces=[.9523 -.3104 0 1.4924 0 0];
                velocidad=wb;
                angulo=.3141;
            case 2
                   V=[0                                                                                       %vector de voltajes inicial
                    0
                    0
                    .01 %.01 original
                    0
                    0];
                enlaces=[1.0002 -.0390 0 8.3698 3.3725 -.0292];
                velocidad=wb;
                angulo=.0014;
                
        end
                
                
        end
        R=[-r  0   0  0  0  0
            0  -r    0  0  0  0
            0  0    -r  0  0  0
            0  0     0   rf  0  0
            0  0     0   0  rs 0
            0  0     0   0   0  rt];
        
        L=[-ld          0           0         mdf   mds     0
            0             -lq        0        0        0         mqt
            0             0           -l0     0        0          0
            -mdf      0           0         lf        lmd          0
            -mds      0           0        lmd      ls         0
            0             -mqt    0         0         0          lt];
        Linv=inv(L);
        variables_estado=[enlaces,velocidad,angulo];
