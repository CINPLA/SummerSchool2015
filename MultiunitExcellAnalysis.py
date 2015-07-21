from __future__ import division
import neo
import numpy as np
import matplotlib.pyplot as plt
import quantities as qt

def plot_waveform(blk, t, u, label):
    unit = get_unit(blk,t,u)
    wave = unit.spiketrains[0].waveforms
    wave = np.mean(wave[:,0,:].T,axis=1)
    wave = normalize(wave,'minmax')
    plt.plot(wave, label=label)
    plt.legend()
    
def normalize(inp, normtype='minmax'):
    if normtype == 'minmax':
        return (inp - inp.min())/np.max(inp - inp.min())
    elif normtype == 'meanstd':
        return (inp - np.mean(inp))/np.std(inp)

def instantaneous_rate(spikes, dt, window, fs=1e4):
    '''
    Here 's' is the spikes, the window size of 'window' is 
    represented by 'dt', window is a function which takes two inputs
    '''
    t = np.linspace(0, spikes.magnitude.max(), fs)*qt.s
    r = np.zeros(t.shape)*qt.Hz
    for s in spikes:
        r += window(dt*qt.s, t - s)
    return r,t

def causalWindow(sig, tau):
    """
    causal window
    """
    a = 1./sig
    w = 0 / tau
    indices = np.where(tau >= 0)
    w[indices] = a**2 * tau[indices] * np.exp(-a * tau[indices]) 
    return w


def gaussianWindow(sig, tau):
    """
    gaussian window
    """
    w = 1.0 / (np.sqrt(2. * np.pi) *  sig) \
    * np.exp(-tau**2 / (2. * sig**2))
    return w
    
    
def rectangularWindow(dt, t):
    """
    rectangular window
    """
    w = 0 / t
    indices = np.where(np.logical_and(-dt/2. <= t, t < dt/2.))
    w[indices] = 1./dt
    return w

def get_unit(blk, t, u):
    for tet in blk.recordingchannelgroups:
        if tet.name == t:
            for unit in tet.units:
                if unit.name == u:
                    return unit

def raster(spikes,ax):
    """
    Raster plot
    """
    for s in spikes.magnitude:
        ax.vlines(s, 0, 1, color = 'b')

    plt.ylim(-.1,1.1)
    ax.set_yticks([])
    ax.set_xticks([])
    return ax