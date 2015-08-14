#!/usr/bin/env python
'''
Cortical pyramidal cell model from Almog and Korngreen (2014) J Neurosci 34:1 182-196
Fast Spiking Basket cell from http://neuromorpho.org/neuroMorpho/neuron_info.jsp?neuron_name=FS-basket

This file is only supplimentary to the Ipython Notebook file with the same name. 
'''

import sys, os
import numpy as np
import pylab as plt
import scipy.fftpack as ff
import neuron
import LFPy

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def exercise1(soma_clamp_params, apic_clamp_params):

    cell_parameters = {
        'morphology': 'A140612.hoc',  # File with cell morphology
        'v_init': -62,
        'passive': False,
        'nsegs_method': None,
        'timeres_NEURON': 2**-3,  # [ms] Should be a power of 2
        'timeres_python': 2**-3,
        'tstartms': -50,  # [ms] Simulation start time
        'tstopms': 50,  # [ms] Simulation end time
        'custom_code': ['cell_model.hoc'] # Loads model specific code
    }

    ### MAKING THE CELL
    cell = LFPy.Cell(**cell_parameters)

    ### MAKING THE INPUT
    apic_clamp_params['idx'] = cell.get_closest_idx(x=-150., y=750., z=0.)
    cl_soma = LFPy.StimIntElectrode(cell, **soma_clamp_params)
    cl_apic = LFPy.StimIntElectrode(cell, **apic_clamp_params)

    cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)

    ### PLOTTING THE RESULTS
    cell_plot_idxs = [soma_clamp_params['idx'], apic_clamp_params['idx']]    
    cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(1./(len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}


    # Plotting the morphology
    plt.figure(figsize=(16,9))
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400], title='Green dots: Inputs')
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    [plt.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=cell_plot_colors[idx], ms=12) for idx in cell_plot_idxs]

    # Plotting the membrane potentials
    plt.subplot(222, title='Membrane potential', xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
    [plt.plot(cell.tvec, cell.vmem[idx, :], c=cell_plot_colors[idx], lw=2) for idx in cell_plot_idxs]

    # Plotting the input currents
    stim_lim = [2*np.min([cl_soma.i, cl_apic.i]), 2*np.max([cl_soma.i, cl_apic.i])]
    plt.subplot(224, title='Input currents', xlabel='Time [ms]', ylabel='nA', ylim=stim_lim)
    plt.plot(cell.tvec, cl_soma.i, c=cell_plot_colors[soma_clamp_params['idx']], lw=2)
    plt.plot(cell.tvec, cl_apic.i, '--', c=cell_plot_colors[apic_clamp_params['idx']], lw=2)

def exercise2(soma_clamp_params, apic_clamp_params, electrode_parameters, noise_level):
    ### MAKING THE CELL
    cell_parameters = {
        'morphology': 'A140612.hoc',
        'v_init': -62,
        'passive': False,
        'nsegs_method': None,
        'timeres_NEURON': 2**-3,  # [ms] Should be a power of 2
        'timeres_python': 2**-3,
        'tstartms': -50,  # [ms] Simulation start time
        'tstopms': 50,  # [ms] Simulation end time
        'custom_code': ['cell_model.hoc'] # Loads model specific code
    }
    cell = LFPy.Cell(**cell_parameters)

    ### MAKING THE INPUT
    apic_clamp_params['idx'] = cell.get_closest_idx(x=-150., y=750., z=0.)
    cl_soma = LFPy.StimIntElectrode(cell, **soma_clamp_params)
    cl_apic = LFPy.StimIntElectrode(cell, **apic_clamp_params)

    cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)

    ### MAKING THE ELECTRODE
    electrode = LFPy.RecExtElectrode(cell, **electrode_parameters)
    electrode.calc_lfp()

    ### PLOTTING THE RESULTS
    cell_plot_idxs = [soma_clamp_params['idx'], apic_clamp_params['idx']]
    cell_idx_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(1./(len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}
    elec_idx_colors = {idx: plt.cm.Reds_r(1./(len(electrode_parameters['x']) + 1) * idx) for idx in range(len(electrode_parameters['x']))}

    plt.figure(figsize=(16,9))
    # Plotting the morphology
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.subplot(141, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400], 
                title='Green dots: Inputs\nRed diamonds: Extracellular electrodes')
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    [plt.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=cell_idx_colors[idx], ms=12) for idx in cell_plot_idxs]
    [plt.plot(electrode_parameters['x'][idx], electrode_parameters['y'][idx], 'D', c=elec_idx_colors[idx], ms=12) for idx in xrange(len(electrode_parameters['x']))]

    # Plotting the membrane potentials
    plt.subplot(242, title='Membrane\npotential', xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
    [plt.plot(cell.tvec, cell.vmem[idx, :], c=cell_idx_colors[idx], lw=2) for idx in cell_plot_idxs]

    # Plotting the input currents
    stim_lim = [2*np.min([cl_soma.i, cl_apic.i]), 2*np.max([cl_soma.i, cl_apic.i])]
    plt.subplot(246, title='Input currents', xlabel='Time [ms]', ylabel='nA', ylim=stim_lim)
    plt.plot(cell.tvec, cl_soma.i, c=cell_idx_colors[soma_clamp_params['idx']], lw=2)
    plt.plot(cell.tvec, cl_apic.i, '--', c=cell_idx_colors[apic_clamp_params['idx']], lw=2)

    # Plotting the extracellular potentials
    LFP = 1000 * electrode.LFP + noise_level * (np.random.random(electrode.LFP.shape) - 0.5)

    plt.subplot(243, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='$\mu$V', xlim=[9, 18])
    [plt.plot(cell.tvec, LFP[idx], c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    norm_LFP = [LFP[idx] - LFP[idx, 0] for idx in xrange(len(electrode_parameters['x']))]
    plt.subplot(247, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='Normalized', xlim=[9, 18])
    [plt.plot(cell.tvec, norm_LFP[idx] / np.max(np.abs(norm_LFP[idx])), c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]


def exercise3(model, input_y_pos):

    ### Making the cell
    if model is 'FS_basket':
        cell_parameters = {
            'morphology': 'FS_basket.hoc',
            'v_init': -65,
            'passive': True,
            'timeres_NEURON': 2**-3,  # Should be a power of 2
            'timeres_python': 2**-3,
            'tstartms': 0,
            'tstopms': 30,
        }
    elif model is 'almog':
        cell_parameters = {
            'morphology': 'A140612.hoc',
            'v_init': -65,
            'passive': True,
            #'nsegs_method': None,
            'timeres_NEURON': 2**-3,  # [ms] Should be a power of 2
            'timeres_python': 2**-3,
            'tstartms': 0,  # [ms] Simulation start time
            'tstopms': 30,  # [ms] Simulation end time
            #'custom_code': ['cell_model.hoc'] # Loads model specific code
        }
    else:
        raise RuntimeError("Wrong model name!")
    cell = LFPy.Cell(**cell_parameters)

    ### Making the synapse
    synapse_params = {
        'idx': cell.get_closest_idx(x=0, y=input_y_pos, z=0),
        'record_current': True,
        'e': 0., #[nA]
        'tau': 2.,
        'weight': 0.005,
        'syntype': 'ExpSyn',
        }
    synapse = LFPy.Synapse(cell, **synapse_params)
    synapse.set_spike_times(np.array([5.]))

    cell.simulate(rec_imem=True, rec_vmem=False)

    ### Finding the time at which to plot the LFP
    time_idx = np.argmax(cell.somav)
    time = cell.tvec[time_idx]

    #  Make dense 2D grid of electrodes
    if model is 'almog':
        x = np.linspace(-500, 500, 30)
        y = np.linspace(-300, 1200, 30)
    elif model is 'FS_basket':
        x = np.linspace(-200, 200, 30)
        y = np.linspace(-200, 200, 30)
    x, y = np.meshgrid(x, y)
    elec_x = x.flatten()
    elec_y = y.flatten()
    elec_z = np.zeros(len(elec_x))
    center_electrode_idx = np.argmin(elec_x**2 + elec_y**2 + elec_z**2)

    electrode_parameters = {
    'sigma': 0.3,              # extracellular conductivity
    'x': elec_x,        # x,y,z-coordinates of contact points
    'y': elec_y,
    'z': elec_z,
    }
    electrode = LFPy.RecExtElectrode(cell, **electrode_parameters)
    electrode.calc_lfp()

    
    plt.figure()
    plt.subplots_adjust(left=0.25)
    
    ### Plotting the cell
    plt.subplot(111, aspect=1, xlabel='x [$\mu$m]', ylabel='y [$\mu$m]', title='Snapshot of LFP at time of maximum\nsomatic membrane potential.\nInput marked by green dot')
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    plt.plot(cell.xmid[cell.synidx], cell.ymid[cell.synidx], 'go', markersize=12)
 
    
    ### Plotting the LFP
    sig_amp = 1000 * electrode.LFP[:, time_idx].reshape(x.shape)
    color_lim = np.max(np.abs(sig_amp))/10
    plt.imshow(sig_amp, origin='lower', extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
               vmin=-color_lim, vmax=color_lim, interpolation='nearest', cmap=plt.cm.bwr_r)
    plt.colorbar(label='$\mu$V')

    plt.axes([0.07, 0.7, 0.1, 0.15], title='Somatic\nmembrane\npotential', xlabel='ms', ylabel='mV',
             xticks=[0, np.max(cell.tvec)], yticks=[cell.somav[0], np.ceil(np.max(cell.somav))], 
             ylim=[cell.somav[0], np.ceil(np.max(cell.somav))])
    plt.plot(cell.tvec, cell.somav)
    plt.plot([cell.tvec[time_idx], cell.tvec[time_idx]], [np.min(cell.somav), np.max(cell.somav)])

def return_freq_and_fft(tvec, sig):
    """ Returns the power and freqency of the input signal"""
    sig = np.array(sig)
    if len(sig.shape) == 1:
        sig = np.array([sig])
    elif len(sig.shape) == 2:
        pass
    else:
        raise RuntimeError("Not compatible with given array shape!")
    timestep = (tvec[1] - tvec[0])/1000.
    sample_freq = ff.fftfreq(sig.shape[1], d=timestep)
    pidxs = np.where(sample_freq >= 0)
    freqs = sample_freq[pidxs]
    Y = ff.fft(sig, axis=1)[:, pidxs[0]]
    amp = np.abs(Y)/Y.shape[1]
    return freqs, amp

def make_white_noise_stimuli(cell, input_idx, max_freq, weight=0.0005):
    """ Makes a white noise input synapse to the cell """ 
    plt.seed(1234)

    # Make an array with sinusoids with equal amplitude but random phases.
    tot_ntsteps = round((cell.tstopms - cell.tstartms) / cell.timeres_NEURON + 1)
    I = np.zeros(tot_ntsteps)
    tvec = np.arange(tot_ntsteps) * cell.timeres_NEURON
    for freq in xrange(1, max_freq + 1):
        I += np.sin(2 * np.pi * freq * tvec/1000. + 2*np.pi*np.random.random())
    input_array = weight * I
    noiseVec = neuron.h.Vector(input_array)

    # Make the synapse
    i = 0
    syn = None
    for sec in cell.allseclist:
        for seg in sec:
            if i == input_idx:
                syn = neuron.h.ISyn(seg.x, sec=sec)
            i += 1
    if syn is None:
        raise RuntimeError("Wrong stimuli index")
    syn.dur = 1E9
    syn.delay = 0 
    noiseVec.play(syn._ref_amp, cell.timeres_NEURON)
    return cell, syn, noiseVec


def exercise4(electrode_parameters, input_y_pos):
    ### MAKING THE CELL
    cell_parameters = {
        'morphology': 'A140612.hoc',
        'v_init': -65,
        'passive': True,
        'timeres_NEURON': 2**-3,  # [ms] Should be a power of 2
        'timeres_python': 2**-3,
        'tstartms': 0,  # [ms] Simulation start time
        'tstopms': 1000,  # [ms] Simulation end time

    }
    cell = LFPy.Cell(**cell_parameters)

    ### MAKING THE INPUT
    max_freq = 500    
    input_idx = cell.get_closest_idx(x=0, y=input_y_pos, z=0.)
    soma_idx = 0
    cell, syn, noiseVec = make_white_noise_stimuli(cell, input_idx, max_freq)

    cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)

    ### MAKING THE ELECTRODE
    electrode = LFPy.RecExtElectrode(cell, **electrode_parameters)
    electrode.calc_lfp()

    ### PLOTTING THE RESULTS

    cell_input_color = 'g'
    elec_idx_colors = {idx: plt.cm.Reds_r(1./(len(electrode_parameters['x']) + 1) * idx) for idx in range(len(electrode_parameters['x']))}

    plt.figure(figsize=(16,9))
    # Plotting the morphology
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.subplot(141, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400], 
                title='Green dots: Inputs\nRed diamonds: Extracellular electrodes')
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    plt.plot(cell.xmid[input_idx], cell.ymid[input_idx], 'o', c=cell_input_color, ms=15)
    [plt.plot(electrode_parameters['x'][idx], electrode_parameters['y'][idx], 'D', c=elec_idx_colors[idx], ms=15) for idx in xrange(len(electrode_parameters['x']))]

    # Plotting the input current
    plt.subplot(242, title='Input current', xlabel='Time [ms]', ylabel='nA')
    plt.plot(cell.tvec, np.array(noiseVec), c=cell_input_color, lw=2)

    freqs, [i_amp] = return_freq_and_fft(cell.tvec, np.array(noiseVec))
    # Plotting the input currents in the frequency domain
    stim_lim = [2*np.min(np.array(noiseVec)), 2*np.max(np.array(noiseVec))]
    plt.subplot(246, title='Input current', xlabel='Frequency [Hz]', ylabel='nA / Hz', xlim=[1, max_freq], ylim=[1e-4, 1e-2])
    plt.loglog(freqs, i_amp, c=cell_input_color, lw=2)

    # Plotting the extracellular potentials
    LFP = 1000 * electrode.LFP
    freqs, LFP_amp = return_freq_and_fft(cell.tvec, LFP)
    plt.subplot(243, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='$\mu$V')
    [plt.plot(cell.tvec, LFP[idx], c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    norm_LFP = [LFP[idx] - LFP[idx, 0] for idx in xrange(len(electrode_parameters['x']))]
    plt.subplot(244, title='Extracellular\npotentials', xlabel='Time [Hz]', ylabel='Normalized')
    [plt.plot(cell.tvec, norm_LFP[idx] / np.max(np.abs(norm_LFP[idx])), c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    # Plotting the extracellular potentials in the frequency domain
    plt.subplot(247, title='Extracellular\npotentials', xlabel='Frequency [Hz]', ylabel='$\mu$V/Hz', xlim=[1, max_freq])
    [plt.loglog(freqs, LFP_amp[idx], c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    plt.subplot(248, title='Extracellular\npotentials', xlabel='Frequency [Hz]', ylabel='Normalized/Hz', xlim=[1, max_freq], ylim=[1e-3, 1e1])
    [plt.loglog(freqs, LFP_amp[idx] / np.max(LFP_amp[idx, 1:]), c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]


