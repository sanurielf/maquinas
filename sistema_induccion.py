# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:30:31 2013

@author: urielsandoval
"""
from copy import copy
from cmath import sin, cos

#import ipdb
import numpy as np
from  numpy.linalg import inv
from numpy import cos, sin, pi
from matplotlib import pyplot as plt
from cmath import sin, cos

from rk4 import rk4
from datos import cargar_datos

class Maquina_Induccion(object):

    def __init__(self, datos, *args, **kwargs):

        self.eventos = []
        self.nombres = {'V': ('V_d', 'V_q', 'V_0', 'V_f'),
                        'I': ('I_d', 'I_q', 'I_0', 'I_f', 'I_s', 'I_t')}
        
        self.Vm = np.ones(6, dtype=float)
        self.V = np.array([0, 1.0, 0, 0, 0, 0], dtype=float)
        self.ws = 2*np.pi*datos['frec']
        
        # Se cargan datos y se crean matrices
        rs = datos['rs']
        rr = datos['rr']

        lr = datos['lr']
        ls = datos['ls']
        lm = datos['lm']

        lss = ls + lm
        lrr = lr + lm

        self.R = np.diag([rs]*3 + [rr]*3)
        
        self.L = np.array([ [lss,    0,   0,  lm,   0, 0], 
                            [0,    lss,   0,   0,  lm, 0],
                            [0,      0,  ls,   0,   0, 0],
                            [lm,     0,   0, lrr,   0, 0],
                            [0,     lm,   0,   0, lrr, 0],
                            [0,      0,   0,   0,   0, lr]], dtype=float)
        self.Linv = inv(self.L)

        self.H = datos['H']
        self.pp = datos['pp']
        self.tm = datos['tm']

        ref = kwargs.get('ref', 'sinc')
        if ref == 'sinc':
            # Referencia síncrona
            self.theta = "0"
            self.wref = "self.ws"
        elif ref == 'est':
            # Referencia estacionaria
            self.theta = "self.ws*t"
            self.wref = "0"
        elif ref == 'rot':
            # Referencia en el rotor
            self.theta = "self.ws*t - x[7]"
            self.wref = "x[6]"

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

            t_eventos[ind].append((True, evento))
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

        Vm = self.Vm
        V = self.V
        R = self.R
        ws = self.ws
        dx = self.dx

        # Calculo de theta de acuerdo a la referenc
        theta = eval(self.theta)
        wref = eval(self.wref)

        # Se calculan corrientes        
        self.I = np.dot(self.Linv, x[:6])
        # Se calcula para eléctrico
        self.te = (x[0] * self.I[1] - x[1] * self.I[0])
        # Inductancias de voltajes rotacionales
        G = np.array([-wref*x[1], wref*x[0], 0, -(wref-x[6])*x[4], (wref-x[6])*x[3], 0])/ws
        # Se calcula el voltaje
        V[0] = Vm[0] * np.sin(theta) 
        V[1] = Vm[1] * np.cos(theta)
        V[2] = 0

        dx[0:6] = (V - np.dot(R, self.I) - G) * ws
        dx[6] = ws *(self.te - self.tm) / (2*self.H)
        dx[7] = x[6]

        return dx
def main():

    global Volt

    dat = cargar_datos('carlos')
    maquina = Maquina_Induccion(dat, ref='sinc')

    enlaces=[0.0]*6
    velocidad = [0]
    angulo = [0]
    maquina.x0 = np.array(enlaces+velocidad+angulo, dtype=float)
    print maquina.x0
    #evento = {'tipo': 'falla_3f',
    #          'ti' : 2,
    #          'tf': 2.01}
    #maquina.eventos.append(evento)          
    maquina.dinamico(5)
    maquina.graficar_todas()
if __name__ == '__main__':
    main()