# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:30:31 2013

@author: urielsandoval
"""
from copy import copy

#import ipdb
import numpy as np
from cmath import sin, pi
from cvxopt import matrix, spdiag


from datos import cargar_datos
from maquina import Maquina
from matrices import inv

class Maquina_Induccion(Maquina):

    def __init__(self, datos, *args, **kwargs):
        super(Maquina_Induccion, self).__init__(datos)

        self.eventos = []
        self.nombres = {'V': ('V_d', 'V_q', 'V_0', 'V_f'),
                        'I': ('I_d', 'I_q', 'I_0', 'I_f', 'I_s', 'I_t')}
        self.V = matrix(0.0, (6, 1))
        self.ws = 2*np.pi*datos['frec']
        
        # Se cargan datos y se crean matrices
        rs = datos['rs']
        rr = datos['rr']

        lr = datos['lr']
        ls = datos['ls']
        lm = datos['lm']

        lss = ls + lm
        lrr = lr + lm

        self.R = spdiag([rs, rs, rs, rr, rr, rr])
        
        self.L = matrix([[lss,    0,   0,  lm,   0, 0], 
                         [0,    lss,   0,   0,  lm, 0],
                         [0,      0,  ls,   0,   0, 0],
                         [lm,     0,   0, lrr,   0, 0],
                         [0,     lm,   0,   0, lrr, 0],
                         [0,      0,   0,   0,   0, lr]]).T
        # Número de enlaces de flujo
        self.nx = self.L.size[0]
        # Se calcula la inversa de L
        self.Linv = inv(self.L)

        self.H = datos['H']
        self.pp = datos['pp']
        self.tm = datos['tm']

        ref = kwargs.get('ref', 'sinc')

        if ref == 'sinc':
            self.theta = "ws*t"
            self.wref =  "1.0"
        elif ref == 'rot':
            self.theta = "thetar"
            self.wref = "wr"
        elif ref == 'est':
            self.theta = "0"
            self.wref = "0"



    def ecuaciones(self, t, x, *args, **kwargs):

        V = self.V
        R = self.R
        ws = self.ws
        dx = self.dx
        # Calculo de theta
        thetar = x[-1]
        wr = x[-2]
        theta = eval(self.theta)
        wref = eval(self.wref)
        # Corrientes de dq0        
        self.I = self.Linv * x[:self.nx]
        # Para eléctrico
        self.te = (x[0] * self.I[1] - x[1] * self.I[0])
        # Inductancias de voltajes rotacionales
        G = matrix([-wref*x[1], wref*x[0], 0, -(wref-wr)*x[4], (wref-wr)*x[3], 0])
        # Voltajes en ABC
        Va = self.VA * sin(self.ws*t).real
        Vb = self.VB * sin(self.ws*t - 2*pi/3).real
        Vc = self.VC * sin(self.ws*t + 2*pi/3).real
        # Voltajes en dq0
        V[:3] = self.transformada_park(theta, matrix([Va, Vb, Vc]))
        # Evaluación de las ED
        dx[0:self.nx] = (V - R * self.I - G) * ws
        dx[-2] = (self.te - self.tm) / (2*self.H)
        dx[-1] = x[6]*ws

        return dx
def main():


    dat = cargar_datos('riaz')
    maquina = Maquina_Induccion(dat, ref='est')

    maquina.x0 = np.zeros(8, dtype=float)

    maquina.eventos.append({'tipo': 'mod_par',
                           'ti': 1.7,
                           'tf': 1.7,
                           'valor': 2})
    # maquina.eventos.append({'tipo': 'falla_3f',
    #                        'ti': 1.7,
    #                        'tf': 1.8,
    #                        'valor': 1})

    maquina.dinamico(5)
    print maquina.X[-1, :]  

    maquina.graficar_todas()
if __name__ == '__main__':
    main()