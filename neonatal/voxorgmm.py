#!/usr/bin/env python

import sys
import numpy as np
from scipy.ndimage.filters import convolve, gaussian_filter
from mrimage import load_mrtrix


def loggauss(X, mu = 0.0, sigma = 1.0, out=None):
    out = np.subtract(X, mu, out=out); out /= sigma   # Z-score
    out = np.square(out, out=out); out *= -0.5        # -0.5 * Z^2
    out -= 0.5 * np.log(2*np.pi) + np.log(sigma)      # offset (normalisation)
    return out


def run_em(X, beta=1.0, sreg=1e-6):
    # initialise
    med = np.median(X)
    mad = np.median(np.abs(X - med)) * 1.4826
    Min, Mout = med, med + 1
    Sin, Sout = mad, mad + 1
    Pin, Pout = 0.9, 0.1
    Win = np.ones_like(X)
    Wout = np.zeros_like(X)
    Rin, Rout = None, None
    # define kernel
    HMRF = np.zeros((3,3,1,1))
    HMRF[0,1] = HMRF[1,0] = HMRF[2,1] = HMRF[1,2] = 1
    # Expectation
    def e_step(X):
        nonlocal Win, Wout, Rin, Rout
        # data energy
        Rin = loggauss(X, Min, Sin, out=Rin); Rin += np.log(Pin)
        Rout = loggauss(X, Mout, Sout, out=Rout); Rout += np.log(Pout)
        # neighbourhood energy
        Rin += beta * convolve(Win, HMRF, mode='wrap')
        Rout += beta * convolve(Wout, HMRF, mode='wrap')
        # calc weights
        lpn = np.logaddexp(Rin, Rout)
        Rin -= lpn; Rout -= lpn;
        Win = np.exp(Rin, out=Win); Win += np.spacing(1)
        Wout = np.exp(Rout, out=Wout); Wout += np.spacing(1)
        return np.mean(lpn);
    # Maximization
    def m_step(X):
        nonlocal Pin, Pout, Min, Mout, Sin, Sout
        Pin = np.mean(Win); Pout = np.mean(Wout)
        Min = np.average(X, weights=Win)
        Mout = np.average(X, weights=Wout)
        Sin = np.sqrt(np.average((X - Min)**2, weights=Win) + sreg)
        Sout = np.sqrt(np.average((X - Mout)**2, weights=Wout) + sreg)
    # Run time !
    ll0 = -np.infty
    for k in range(50):
        ll = e_step(X)
        m_step(X)
        if (np.abs(ll - ll0) < 1e-3):
            break
        ll0 = ll
    # Output
    print("mean: {:.4f}, {:.4f} | stdev: {:.4f}, {:.4f}".format(Min, Mout, Sin, Sout))
    return Win


if __name__ == '__main__':
    I = load_mrtrix(sys.argv[1])
    X = I.data.astype(np.float32)
    I.data = run_em(X)
    I.save(sys.argv[2])


