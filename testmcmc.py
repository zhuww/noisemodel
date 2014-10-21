import numpy as np
import scipy.linalg as sl
import numpy.random as nr
from pylab import *
from ProgressBar import progressBar as PB
#mvn = np.random.multivariate_normal
#print cov
#cov_I = np.inv(cov)


def Metropolis(func, x0, m=100, n=1000, ss = 1.0):
    f0 = func(x0)
    fmin = f0
    xmin = x0
    x = x0
    result = []
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
            result.append([f0] + list(x))

    return result, xmin, fmin


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
    


def SliceSampleMC(func, x0, m=100, n=1000, ordered=True, ss=1.0, progressbar=False):
    lpar = len(x0)
    w = np.abs(x0) * ss #
    w[w==0.]=1.#prevent w from being 0
    x = np.copy(x0)
    f0 = func(x0)
    Pn = lambda x:exp(func(x) - f0) 
    #normalize the loglikihood to f0
    p = 1.
    xmax = np.copy(x0)
    pmax = p
    result = np.array(np.hstack((p, x0)))
    counter = 0
    if progressbar:pb = PB(maxValue = m+n, totalWidth=45)
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
            if p > pmax:
                xmax = np.copy(x)
                pmax = np.float(p)
            if counter > m:
                result = np.vstack((result, np.hstack((p,x))))
            if progressbar and (counter % 10) == 0:pb(counter)
            counter += 1
            #print counter, j+1
    return result[1:,...], xmax, pmax


if __name__ == "__main__":
    from scipy.interpolate import griddata
    m = 3
    N = np.zeros((m,m))
    np.fill_diagonal(N, np.ones(m))
    cov = 0.3 * nr.randn(m, m) + N
    cf = sl.cho_factor(cov)
    def func(p):
        return -0.5 * np.dot(p.T, sl.cho_solve(cf, p)) - 0.5 * 2.*np.sum(np.log(np.diag(cf[0])))

    x0 = nr.randn(m)
    #res, xmax, pmax = SliceSampleMC(func, x0, n = 10000, ss=0.9, ordered=False) 
    res, xmax, pmax = Metropolis(func, x0, n=10000, ss=1.) 
    res = np.array(res)
    parray = res[:,0]
    res[:,0] = (parray - pmax)
    print res.shape
    x = res[:,1]
    y = res[:,2]
    z = res[:,0]
    xi = np.linspace(np.min(x), np.max(x), 30)
    yi = np.linspace(np.min(y), np.max(y), 30)
    zi= griddata((x,y), z, (xi[None,:], yi[:,None]), method='cubic')
    #plot(z)
    print np.max(z), np.min(z)
    levels = [0.5*np.min(z)]
    CS = contour(xi,yi,zi, levels)
    imshow(zi, cmap=cm.gray, extent=(np.max(x), np.min(x), np.min(y),np.max(y)))
    #scatter(x, y)
    #contour(res[:,1], res[:2], res[:0])
    show()
