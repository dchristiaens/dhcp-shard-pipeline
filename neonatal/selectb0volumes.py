#!/usr/bin/env python

import sys
import numpy as np
from mrimage import load_mrtrix
from skimage.restoration import unwrap_phase
from scipy.ndimage import convolve1d


def filterphasemap(image):
    uphase = np.zeros(image.shape)
    for k in range(image.shape[-1]):
        p = np.squeeze(image[:,:,:,k])
        uphase[:,:,:,k] = unwrap_phase(p)
    return convolve1d(uphase, [-1,2,-1], axis=2)


def calcscore(image, mask):
    return np.mean(np.abs(image[mask]), axis=0)


if __name__ == '__main__':
    I = load_mrtrix(sys.argv[1])
    mask = np.squeeze(load_mrtrix(sys.argv[2]).data > 0.5)
    petable = np.mod(np.arange(I.shape[3]), 4)  # hard code PE scheme as 4-axis rotation
    b0idx = np.squeeze(np.argwhere(I.grad[:,3]<1))
    b0pe = petable[b0idx]
    fphase = filterphasemap(I.data[:,:,:,b0idx])
    score = calcscore(fphase, mask)
    argidx = np.zeros((8,), dtype=int)
    for k, pe in enumerate(range(4)):
        j = np.argsort(score[b0pe==pe])[:2]
        argidx[2*k:2*k+2] = np.squeeze(np.argwhere(b0pe==pe))[j]
    print(','.join([str(i) for i in b0idx[argidx]]))



