import numpy as np
import scipy.linalg as sl
import numpy.random as nr
from pylab import *
#mvn = np.random.multivariate_normal
#print cov
#cov_I = np.inv(cov)


def Metropolis(func, x0, m=100, n=1000, ss = 1.0):
    result = []
    f0 = func(x0)
    fmin = f0
    xmin = x0
    x = x0
    plen = len(x0)
    for counter in range(m+n):
        xn = x + nr.randn(plen) * np.abs(x0) * ss
        f = func(xn)
        if f > fmin:
            fmin = f
            xmin = xn
        df = f - f0
        unif = nr.uniform(0,1,1)
        if unif < exp(df):
            x = xn
            f0 = f
        if counter >= m:
            result.append(list(x))

    return result


def SliceSampler(Px, x0, u, w):
    xmax = x0+w
    xmin = x0-w
    while Px(xmax) > u:
        xmax += w
    while Px(xmin) > u:
        xmin -= w
    x = nr.uniform(xmin, xmax)
    p = Px(x)
    while p < u:
        #print Px(x), u, w, x, x0, xmin, xmax
        if x > x0:
            xmax = x
        else:
            xmin = x
        x = nr.uniform(xmin, xmax)
        p = Px(x)
    #return x, p, max((xmax-xmin)/2, 0.01)
    return x, p, (xmax-xmin)/2
    


def SliceSampleMC(func, x0, m=100, n=1000, ordered=True, ss=1.0):
    result = np.array(x0)
    lpar = len(x0)
    w = np.abs(x0) * ss #
    w[w==0.]=1.#prevent w from being 0
    x = np.copy(x0)
    f0 = func(x0)
    Pn = lambda x:exp(func(x) - f0) 
    #normalize the loglikihood to f0
    p = 1.
    counter = 0
    while counter < (m+n):
        for i in range(lpar):
            if ordered:
                j = i
            else:
                j = nr.choice(lpar)
            u = nr.uniform(0, p)
            def Px(z, x=x, j=j):
                y = np.copy(x) 
                y[j] = z
                return Pn(y)
            xn,p,wn = SliceSampler(Px, x[j], u, w[j])
            x[j] = xn
            w[j] = wn
            #print x[j], xn
            if counter > m:
                result = np.vstack((result, x))
            counter += 1
            print counter, j+1
    return result[1:,...]


if __name__ == "__main__":
    m = 3
    N = np.zeros((m,m))
    np.fill_diagonal(N, np.ones(m))
    cov = 0.3 * nr.randn(m, m) + N
    cf = sl.cho_factor(cov)
    def func(p):
        return -0.5 * np.dot(p.T, sl.cho_solve(cf, p)) - 0.5 * 2.*np.sum(np.log(np.diag(cf[0])))

    x0 = nr.randn(m)
    res = SliceSampleMC(func, x0, ss=0.1) 
    #res = Metropolis(func, x0)#, n=10000)#, ss=0.8) 
    res = np.array(res)
    print res.shape
    scatter(res[:,1], res[:,2])
    show()
