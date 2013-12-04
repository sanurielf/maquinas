import numpy as np


from sistema_sincrona import Maquina_Sincrona
from datos import cargar_datos


def krauze(simulacion = 1):
    if simulacion ==1 :
        dat = cargar_datos('krauze_hidro')
        cond_ini = [0.9999991096401809, 1.7487637739538667e-08, 0.0, 1.533434892167535, 1.084575870197835, 1.3022709749357489e-08, 376.9911159870241, -0.001339224497901841]
        
        maquina = Maquina_Sincrona(dat)
        maquina.x0 = np.array(cond_ini)
        
        maquina.eventos.append({'tipo':'mod_par','valor': 0.5, 'ti': 0, 'tf':100})
        
        maquina.dinamico(15)
        maquina.graficar_todas()
    elif simulacion ==2 :
        dat = cargar_datos('krauze_hidro')
        cond_ini = [0.9526469972241254, -0.3093381597437568, 0.0, 1.4936821504967202, 1.0442214825114955, -0.23200361964778743, 376.9911184367233, 0.31297264115349566]

        
        maquina = Maquina_Sincrona(dat)
        maquina.x0 = np.array(cond_ini)
        
        maquina.eventos.append({'tipo':'falla_3f', 'ti': .1, 'tf':.56})
        
        maquina.dinamico(1)
        #maquina.graficar_todas()

def main():
    krauze(2)


if __name__ == '__main__':
    main()