# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:30:31 2013

@author: urielsandoval
"""
from copy import copy


import numpy as np
from  numpy.linalg import inv
from numpy import cos, sin, pi, sqrt
from matplotlib import pyplot as plt


from rk4 import rk4
from datos import cargar_datos

Volt = []

T1 = 2./3 * np.array([[1, -0.5, -0.5],
                      [0, sqrt(3)/2, -sqrt(3)/2]])

defase = 2*pi/3
raiz1_2 = 1/sqrt(2)
raiz2_3 = sqrt(2./3)

class Maquina(object):

    def __init__(self, datos):

        # self.eventos = []
        # self.nombres = {'V': ('V_d', 'V_q', 'V_0', 'V_f'),
        #                 'I': ('I_d', 'I_q', 'I_0', 'I_f', 'I_s', 'I_t')}
        # self.Vm = np.ones(6)
        # self.V = np.array([0, 0, 0, datos['Vfd'], 0, 0])
        # self.ws = 2*np.pi*datos['frec']
        
        # # Se cargan datos y se crean matrices
        # r = datos['r']
        # rf = datos['rf']
        # rs = datos['rs']
        # rt = datos['rt']
        # ld = datos['ld']
        # lq = datos['lq']
        # l0 = datos['l0']
        # ll = datos['ll']
        # llfd = datos['llfd']
        # llkd = datos['llkd']
        # llkq = datos['llkq']
        # lmd = ld - ll
        # lmq = lq - ll
        # lf = llfd + lmd
        # ls = llkd + lmd
        # lt = llkq + lmq

        # mdf = mds = lmd
        # mqt = lmq

        # self.R = np.diag([-r, -r, -r, rf, rs, rt])
        # self.L = np.array([[-ld, 0, 0, mdf, mds, 0], 
        #                     [0, -lq, 0, 0, 0, mqt],
        #                     [0, 0, -l0, 0, 0, 0],
        #                     [-mdf, 0, 0, lf, lmd, 0],
        #                     [-mds, 0, 0, lmd, ls, 0],
        #                     [0, -mqt, 0, 0, 0, lt]])
        # self.Linv = inv(self.L)

        # self.H = datos['H']
        # self.pp = datos['pp']
        # self.tm = datos['tm']
        self.VA = 1
        self.VB = 1
        self.VC = 1

    def abc2ab(self, Vabc):

        return np.dot(T1, Vabc)

    def ab2dq(self, beta, Vab):

        T = np.array([[cos(beta), sin(beta)],
                      [-sin(beta), cos(beta)]])

        return np.dot(T, Vab)

    def transformada_park(self, theta, Vabc):

        C = raiz2_3 * np.array([[cos(theta), cos(theta - defase), cos(theta + defase)],\
                            [-sin(theta), -sin(theta - defase), -sin(theta + defase)],\
                            [raiz1_2, raiz1_2, raiz1_2]])


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
        self.data['w'][0] = self.x0[6] #* 120 / (2*np.pi*self.pp)
        self.data['d'][0] = self.x0[7] #* 180 / np.pi


        for p in xrange(1, N):
            #print 'Iteracion %d'%p
            self.p = p
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

        #self.data['w'] *= 120 / (2*np.pi*self.pp)
        #self.data['d'] *= 180 / np.pi

        #print self.X[-1,:].tolist()


    def graficar_todas(self):
        plt.close('all')
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
        ax1.plot(self.t, self.data['w'], linewidth=1, label='$\omega$')
        ax1.legend(loc=0)
        ax2.plot(self.t, self.data['d'], linewidth=1, label='$\delta$')
        ax2.legend(loc=0)
        [ax3.plot(self.t, self.data['V'][:, p], linewidth=1, label='$'+nombre+'$') for p,nombre in zip(range(4),self.nombres['V'])]
        ax3.legend(loc=0)
        [ax4.plot(self.t, self.data['I'][:, p], linewidth=1, label='$'+nombre+'$') for p,nombre in zip(range(6), self.nombres['I'])]
        ax4.legend(loc=0)

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

