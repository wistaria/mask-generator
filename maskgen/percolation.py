# coding:utf-8

if __name__ == 'maskgen.percolation':
    from . import engine
else:
    import engine

def percolation(seed, L, n, prob = -1.0):
    sim = engine.percolation(seed, L, prob)
    for i in range(n):
        sim.update()
        yield sim.config()
