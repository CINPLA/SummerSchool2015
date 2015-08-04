#!/usr/bin/env python
'''
Multi-compartment NEURON model from Almog and Korngreen (2014) J Neurosci 34:1 182-196
'''

import numpy as np
import pylab as plt
import LFPy
import sys, os

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def exercise1(cell_parameters, soma_clamp_params, apic_clamp_params):
    ### MAKING THE CELL
    cell = LFPy.Cell(**cell_parameters)

    ### MAKING THE INPUT
    apic_clamp_params['idx'] = cell.get_closest_idx(x=-150., y=750., z=0.)

    cl_soma = LFPy.StimIntElectrode(cell, **soma_clamp_params)
    cl_apic = LFPy.StimIntElectrode(cell, **apic_clamp_params)

    cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)

    ### PLOTTING THE RESULTS
    cell_plot_idxs = [soma_clamp_params['idx'], apic_clamp_params['idx']]
    cell_plot_colors = {soma_clamp_params['idx']: 'b',
                        apic_clamp_params['idx']: 'g'}

    # Plotting the morphology
    plt.figure(figsize=(16,9))
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400])
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    [plt.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=cell_plot_colors[idx], ms=7) for idx in cell_plot_idxs]

    # Plotting the membrane potentials
    plt.subplot(222, title='Membrane potential', xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
    [plt.plot(cell.tvec, cell.vmem[idx, :], c=cell_plot_colors[idx], lw=2) for idx in cell_plot_idxs]

    # Plotting the input currents
    stim_lim = [2*np.min([cl_soma.i, cl_apic.i]), 2*np.max([cl_soma.i, cl_apic.i])]
    plt.subplot(224, title='Input currents', xlabel='Time [ms]', ylabel='nA', ylim=stim_lim)
    plt.plot(cell.tvec, cl_soma.i, c=cell_plot_colors[soma_clamp_params['idx']], lw=2)
    plt.plot(cell.tvec, cl_apic.i, '--', c=cell_plot_colors[apic_clamp_params['idx']], lw=2)

def exercise2(cell_parameters, soma_clamp_params, apic_clamp_params, electrode_parameters):
    ### MAKING THE CELL
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
    cell_idx_colors = {cell_plot_idxs[idx]: plt.cm.Blues_r(1./(len(cell_plot_idxs) + 1) * idx) for idx in range(len(cell_plot_idxs))}
    elec_idx_colors = {idx: plt.cm.Reds_r(1./(len(electrode_parameters['x']) + 1) * idx) for idx in range(len(electrode_parameters['x']))}

    plt.figure(figsize=(16,9))
    # Plotting the morphology
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.subplot(141, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400])
    [plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
    [plt.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=cell_idx_colors[idx], ms=7) for idx in cell_plot_idxs]
    [plt.plot(electrode_parameters['x'][idx], electrode_parameters['y'][idx], 'D', c=elec_idx_colors[idx], ms=7) for idx in xrange(len(electrode_parameters['x']))]

    # Plotting the membrane potentials
    plt.subplot(242, title='Membrane\npotential', xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
    [plt.plot(cell.tvec, cell.vmem[idx, :], c=cell_idx_colors[idx], lw=2) for idx in cell_plot_idxs]

    # Plotting the input currents
    stim_lim = [2*np.min([cl_soma.i, cl_apic.i]), 2*np.max([cl_soma.i, cl_apic.i])]
    plt.subplot(246, title='Input currents', xlabel='Time [ms]', ylabel='nA', ylim=stim_lim)
    plt.plot(cell.tvec, cl_soma.i, c=cell_idx_colors[soma_clamp_params['idx']], lw=2)
    plt.plot(cell.tvec, cl_apic.i, '--', c=cell_idx_colors[apic_clamp_params['idx']], lw=2)

    # Plotting the extracellular potentials
    plt.subplot(243, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='mV', xlim=[9, 18])
    [plt.plot(cell.tvec, electrode.LFP[idx], c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    norm_LFP = [electrode.LFP[idx] - electrode.LFP[idx, 0] for idx in xrange(len(electrode_parameters['x']))]
    plt.subplot(247, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='Normalized', xlim=[9, 18])
    [plt.plot(cell.tvec, norm_LFP[idx] / np.max(np.abs(norm_LFP[idx])), c=elec_idx_colors[idx], lw=2) for idx in xrange(len(electrode_parameters['x']))]

    LFP_amplitude_decay = np.max(electrode.LFP, axis=1) - np.min(electrode.LFP, axis=1)

    plt.subplot(244, title='Extracellular\npotential decay', xlabel='x [$\mu m$]', ylabel='mV')
    plt.plot(-electrode_parameters['x'], LFP_amplitude_decay)
    plt.grid(True)

    plt.subplot(248, title='Extracellular\npotential decay', xlabel='x [$\mu m$]', ylabel='Normalized', ylim=[0, 1.1])
    plt.plot(-electrode_parameters['x'], LFP_amplitude_decay / np.max(LFP_amplitude_decay))
    plt.grid(True)
