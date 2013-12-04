# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:30:31 2013

@author: urielsandoval
"""
from copy import copy


import numpy as np
from  numpy.linalg import inv
from numpy import cos, sin, pi
from matplotlib import pyplot as plt


from rk4 import rk4
from datos import cargar_datos

Volt = []
class Maquina_Sincrona(object):

    def __init__(self, datos):

        self.eventos = []
        self.nombres = {'V': ('V_d', 'V_q', 'V_0', 'V_f'),
                        'I': ('I_d', 'I_q', 'I_0', 'I_f', 'I_s', 'I_t')}
        self.Vm = np.ones(6)
        self.V = np.array([0, 0, 0, datos['Vfd'], 0, 0])
        self.ws = 2*np.pi*datos['frec']
        
        # Se cargan datos y se crean matrices
        r = datos['r']
        rf = datos['rf']
        rs = datos['rs']
        rt = datos['rt']
        ld = datos['ld']
        lq = datos['lq']
        l0 = datos['l0']
        ll = datos['ll']
        llfd = datos['llfd']
        llkd = datos['llkd']
        llkq = datos['llkq']
        lmd = ld - ll
        lmq = lq - ll
        lf = llfd + lmd
        ls = llkd + lmd
        lt = llkq + lmq

        mdf = mds = lmd
        mqt = lmq

        self.R = np.diag([-r, -r, -r, rf, rs, rt])
        self.L = np.array([[-ld, 0, 0, mdf, mds, 0], 
                            [0, -lq, 0, 0, 0, mqt],
                            [0, 0, -l0, 0, 0, 0],
                            [-mdf, 0, 0, lf, lmd, 0],
                            [-mds, 0, 0, lmd, ls, 0],
                            [0, -mqt, 0, 0, 0, lt]])
        self.Linv = inv(self.L)

        self.H = datos['H']
        self.pp = datos['pp']
        self.tm = datos['tm']
        self.VA = 1
        self.VB = 1
        self.VC = 1

    def calcula_voltajes(self, theta, Vabc):



        C = 2./3 * np.array([[cos(theta), cos(theta - 2*pi/3), cos(theta + 2*pi/3)],\
                            [-sin(theta), -sin(theta - 2*pi/3), -sin(theta + 2*pi/3)],\
                            [0.5, 0.5, 0.5]])


        print C
        print Vabc
        print 'multiplicacion'
        print np.dot(C, Vabc)
        return np.dot(C, Vabc)

    def aplicar_evento(self, aplicar,evento):

        if aplicar:
            print 'Se apica evento'
            if evento['tipo'] == 'falla_3f':
                self.Vm1 = copy(self.Vm)
                self.VA = 0
                self.VB = 0
                self.VC = 0
                self.Vm[0] = 0
                self.Vm[1] = 0
            elif evento['tipo'] == 'mod_par':
                self.tm = evento['valor']

        else:
            print 'Se libera evento'
            if evento['tipo'] == 'falla_3f':
                self.VA = 1
                self.VB = 1
                self.VC = 1
                self.Vm[0] = self.Vm1[0]
                self.Vm[1] = self.Vm1[0]
            elif evento['tipo'] == 'mod_par':
                pass

    def crea_tiempo_eventos(self, t):

        t_eventos = [[] for _ in range(t.size)]
        for evento in self.eventos:
            ind = (np.abs(t - evento['ti'])).argmin()
            ind2 = (np.abs(t - evento['tf'])).argmin()

            t_eventos[ind].append((True,evento))
            t_eventos[ind2].append((False, evento))

        return t_eventos

    def dinamico(self, tf, **kwargs):

        h = kwargs.get('h', 0.001)
        N = int((tf) / h)
        self.t = np.linspace(0, tf, N)
        ceros = np.zeros(N)
        ceros_n = np.zeros((N, 8))
        ceros_n_2 = np.zeros((N, 6))


        t_eventos = self.crea_tiempo_eventos(self.t)
        self.X = ceros_n.copy()

        self.data = {'V': ceros_n_2.copy(),
                    'I': ceros_n_2.copy(),
                    'Te': ceros.copy(),
                    'Tm': ceros.copy(),
                    'w': ceros.copy(),
                    'd': ceros.copy()}

        self.dx = np.zeros(8)

        self.X[0, :] = self.x0
        self.data['w'][0] = self.x0[6] * 120 / (2*np.pi*self.pp)
        self.data['d'][0] = self.x0[7] * 180 / np.pi


        for p in xrange(1, N):
            #print 'Iteracion %d'%p
            if t_eventos[p-1]:
                for evento in t_eventos[p-1]:
                    print evento
                    self.aplicar_evento(evento[0], evento[1])



            self.X[p, :] = rk4(self.ecuaciones, self.X[p-1, :], self.t[p], h)
            # Se almacenan variables
            self.data['V'][p] = self.V
            self.data['w'][p] = self.X[p, 6]
            self.data['d'][p] = self.X[p, 7]
            self.data['I'][p] = self.I
            self.data['Te'][p] = self.te
            self.data['Tm'][p] = self.tm

        self.data['w'] *= 120 / (2*np.pi*self.pp)
        self.data['d'] *= 180 / np.pi

        print self.X[-1,:].tolist()


    def graficar_todas(self):
        plt.close('all')
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
        ax1.plot(self.t, self.data['w'], linewidth=1, label='$\omega$',)
        ax1.legend(loc=1)
        ax2.plot(self.t, self.data['d'], linewidth=1, label='$\delta$')
        ax2.legend(loc=1)
        [ax3.plot(self.t, self.data['V'][:, p], linewidth=1, label='$'+nombre+'$') for p,nombre in zip(range(4),self.nombres['V'])]
        ax3.legend(loc=1)
        [ax4.plot(self.t, self.data['I'][:, p], linewidth=1, label='$'+nombre+'$') for p,nombre in zip(range(6), self.nombres['I'])]
        ax4.legend(loc=1)

        f1 = plt.figure()

        ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
        ax2 = plt.subplot2grid((2, 2), (0, 0))
        ax1 = plt.subplot2grid((2, 2), (0, 1))

        ax1.plot(self.data['d'], self.data['Te'], linewidth=1)
        ax1.set_xlabel('$\delta$')
        ax1.set_ylabel('$T_e$')

        ax2.plot(self.data['w'], self.data['d'], linewidth=1)
        ax2.set_xlabel('$\omega_r$')
        ax2.set_ylabel('$\delta$')

        ax3.plot(self.t, self.data['Te'], linewidth=1, label='$T_e$')
        ax3.plot(self.t, self.data['Tm'], linewidth=1, label='$T_m$')
        ax3.set_xlabel('$t$')
        ax3.set_ylabel('Par')

        ax3.legend(loc=1)

        plt.show()

    def ecuaciones(self, t, x):
        global Volt

        Vm = self.Vm
        V = self.V
        R = self.R
        ws = self.ws
        dx = self.dx

        # Se calculan corrientes        
        self.I = np.dot(self.Linv, x[:6])
        # Se calcula para eléctrico
        self.te = x[0] * self.I[1] - x[1] * self.I[0]
        # Inductancias de voltajes rotacionales
        G = np.array([x[1], -x[0], 0, 0, 0, 0])
        # Se calculan voltaje
        alfa = pi/2 - x[7]

        Vabc = np.array([self.VA* sin(0*ws*t + alfa), self.VB * sin(0*ws*t - 2*pi/3 + alfa), self.VC * sin(0*ws*t + 2*pi/3 + alfa)])
        Vdq0 = self.calcula_voltajes(0, Vabc)
        Volt.append(Vdq0)
        #print Vdq0[0] 
        V[0] = Vm[0] * np.sin(x[7]) #Vdq0[0] #
        #print V[0]
        V[1] = Vm[1] * np.cos(x[7]) #Vdq0[1] #
        V[2] = Vdq0[2]

        dx[0:6] = (V - np.dot(R, self.I) + (x[6]/ws) * G) * ws
        dx[6] = ws *(self.tm - self.te) / (2*self.H)
        dx[7] = x[6] - ws

        return dx
def main():

    global Volt

    dat = cargar_datos('kundur')
    maquina = Maquina_Sincrona(dat)

    enlaces=[.9523, -.3104, 0, 1.4924, 0, 0]
    #enlaces = [1.0002, -.0390, 0, 8.3698, 3.3725, -.0292]
    velocidad = [2*np.pi*60]
    angulo = [0.3141]
    #angulo = [.0014]    
    maquina.x0 = np.array(enlaces+velocidad+angulo)

    evento = {'tipo': 'falla_3f',
              'ti' : 2,
              'tf': 2.01}
    maquina.eventos.append(evento)          
    maquina.dinamico(0.016)
    #maquina.graficar_todas()
    plt.plot(Volt)
    plt.show()
if __name__ == '__main__':
    main()