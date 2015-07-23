#!/usr/bin/env python
'''
Multi-compartment NEURON model from Almog and Korngreen (2014) J Neurosci 34:1 182-196
'''

import numpy as np
import pylab as plt
import LFPy

### MAKING THE CELL
cell_parameters = {
    'morphology': 'A140612.hoc',
    'v_init': -62,
    'passive': False,
    'nsegs_method': None,
    'timeres_NEURON': 2**-3,  # Should be a power of 2
    'timeres_python': 2**-3,
    'tstartms': -50,
    'tstopms': 50,
    'custom_code': ['cell_model.hoc']
}
cell = LFPy.Cell(**cell_parameters)

### MAKING THE INPUT
soma_idx = 0
apic_idx = cell.get_closest_idx(x=-150., y=750., z=0.)

soma_clamp_params = {
    'idx': soma_idx,
    'record_current': True,
    'amp': -0.5, #[nA]
    'dur': 30.,
    'delay': 10,
    'pptype': 'IClamp',
    }

apic_clamp_params = {
    'idx': apic_idx,
    'record_current': True,
    'amp': -.25, #[nA]
    'dur': 30.,
    'delay': 10,
    'pptype': 'IClamp',
    }
cl_soma = LFPy.StimIntElectrode(cell, **soma_clamp_params)
cl_apic = LFPy.StimIntElectrode(cell, **apic_clamp_params)

cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)

### PLOTTING THE RESULTS
cell_plot_idxs = [soma_idx, apic_idx]
cell_plot_colors = {soma_idx: 'b',
                    apic_idx: 'g'}

# Plotting the morphology
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
plt.plot(cell.tvec, cl_soma.i, c=cell_plot_colors[soma_idx], lw=2)
plt.plot(cell.tvec, cl_apic.i, '--', c=cell_plot_colors[apic_idx], lw=2)

plt.savefig('exercise_1.png')

