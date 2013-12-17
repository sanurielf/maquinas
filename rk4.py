
def rk4(funcion, x, t, h):

    f1 = h * funcion(t, x)
    #print x.size
    #print f1.size
    f2 = h * funcion(t + h / 2, x + f1 / 2)
    f3 = h * funcion(t + h / 2, x + f2 / 2)
    f4 = h * funcion(t + h, x + f3)
    x1 = x + (f1 + 2 * (f2 + f3) + f4) / 6
    return x1
