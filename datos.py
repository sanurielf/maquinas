
def cargar_datos(caso):

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
        tm = 0.85
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