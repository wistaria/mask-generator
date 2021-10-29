# coding:utf-8

if __name__ == 'maskgen.ising':
    from . import engine
else:
    import engine

def ising(seed, L, n, beta = -1.0, therm = 128, interval = 128):
    sim = engine.ising(seed, L, beta)
    for i in range(therm):
        sim.update()
    for i in range(n):
        if (i > 0):
            for j in range(interval + 1):
                sim.update()
        yield sim.config()
