
def cargar_datos(caso):

    if caso == 'carlos':

        datos = {}
        datos['rs'] = 0.01965
        datos['rr'] = 0.01909
        datos['ls'] = 0.0397
        datos['lr'] = 0.0397
        datos['lm'] = 1.354
        datos['H'] =  0.09526*2
        datos['pp'] = 2
        datos['frec'] = 60
        datos['tm'] = 1
        return datos
    
    if caso == 'krause_motor':

        datos = {}
        datos['rs'] = 0.0453
        datos['rr'] = 0.0222
        datos['ls'] = 0.0775
        datos['lr'] = 0.0322
        datos['lm'] = 2.042
        datos['H'] =  0.5
        datos['pp'] = 3
        datos['frec'] = 60
        datos['tm'] = 0.8
        return datos

    if caso=='krauze_hidro':

        # Parametros resistencias
        r = 0.0019
        rf = 4.1e-4
        rs = 1.41e-2
        rt = 1.36e-2
        # Parametros inductancias
        ll = 0.12
        llfd = .2049
        llkd = .16
        llkq = .1029

        ld = .85
        lq = .48
        l0 = .12

        H = 7.5
        pp = 64
        tm = 0.0
        Vfd = 0.00089936
        frec = 60

    elif caso=='kundur':
        r = 0.003
        rf = 0.0006
        rs = 0.0284
        rt = 0.00619

        ll = 0.15
        llfd = 0.165
        llkd = 0.1713
        llkq = 0.7252

        lad = 1.66
        laq = 1.61
        ld = lad + ll
        lq = laq + ll
        l0 = ld/2

        H = 2.5
        pp = 1
        tm = 0.85
        Vfd = 0.00089936
        frec = 60

    elif caso=='krauze':


        pass


    datos = {'r': r,
            'rf': rf,
            'rs': rs,
            'rt': rt,
            'll': ll,
            'llfd': llfd,
            'llkd': llkd,
            'llkq': llkq,
            'ld': ld,
            'lq': lq,
            'l0': l0,
            'H': H,
            'pp': pp,
            'tm': tm,
            'Vfd': Vfd,
            'frec': frec}
    return datos