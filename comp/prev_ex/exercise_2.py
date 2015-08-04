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
    'timeres_NEURON': 2**-4,  # Should be a power of 2
    'timeres_python': 2**-4,
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
    'amp': .5, #[nA]
    'dur': 20.,
    'delay': 10,
    'pptype': 'IClamp',
    }

apic_clamp_params = {
    'idx': apic_idx,
    'record_current': True,
    'amp': 0.5, #[nA]
    'dur': 20.,
    'delay': 10,
    'pptype': 'IClamp',
    }
cl_soma = LFPy.StimIntElectrode(cell, **soma_clamp_params)
cl_apic = LFPy.StimIntElectrode(cell, **apic_clamp_params)

cell.simulate(rec_imem=True, rec_vmem=True, rec_istim=True)


### MAKING THE ELECTRODE
elec_x = -25. * np.ones(3)
elec_y = np.array([0., 500., 800.])
elec_z = np.zeros(len(elec_x))
electrode_parameters = {
    'sigma': 0.3,              # extracellular conductivity
    'x': elec_x,        # x,y,z-coordinates of contact points
    'y': elec_y,
    'z': elec_z,
}
electrode = LFPy.RecExtElectrode(cell, **electrode_parameters)
electrode.calc_lfp()




### PLOTTING THE RESULTS
cell_plot_idxs = [soma_idx, apic_idx]
cell_idx_colors = {cell_plot_idxs[idx]: plt.cm.Blues_r(1./(len(cell_plot_idxs) + 1) * idx) for idx in range(len(cell_plot_idxs))}
elec_idx_colors = {idx: plt.cm.Reds_r(1./(len(elec_x) + 1) * idx) for idx in range(len(elec_x))}

plt.figure(figsize=[12, 6])
# Plotting the morphology
plt.subplots_adjust(hspace=0.5, wspace=0.5)
plt.subplot(141, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]', xlim=[-400, 400], xticks=[-400, 0, 400])
[plt.plot([cell.xstart[idx], cell.xend[idx]], [cell.ystart[idx], cell.yend[idx]], c='k') for idx in xrange(cell.totnsegs)]
[plt.plot(cell.xmid[idx], cell.ymid[idx], 'o', c=cell_idx_colors[idx], ms=7) for idx in cell_plot_idxs]
[plt.plot(elec_x[idx], elec_y[idx], 'D', c=elec_idx_colors[idx], ms=7) for idx in xrange(len(elec_x))]

# Plotting the membrane potentials
plt.subplot(242, title='Membrane\npotential', xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
[plt.plot(cell.tvec, cell.vmem[idx, :], c=cell_idx_colors[idx], lw=2) for idx in cell_plot_idxs]

# Plotting the input currents
stim_lim = [2*np.min([cl_soma.i, cl_apic.i]), 2*np.max([cl_soma.i, cl_apic.i])]
plt.subplot(246, title='Input currents', xlabel='Time [ms]', ylabel='nA', ylim=stim_lim)
plt.plot(cell.tvec, cl_soma.i, c=cell_idx_colors[soma_idx], lw=2)
plt.plot(cell.tvec, cl_apic.i, '--', c=cell_idx_colors[apic_idx], lw=2)

# Plotting the extracellular potentials
plt.subplot(243, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='mV', xlim=[9, 18])
[plt.plot(cell.tvec, electrode.LFP[idx], c=elec_idx_colors[idx], lw=2) for idx in xrange(len(elec_x))]

norm_LFP = [electrode.LFP[idx] - electrode.LFP[idx, 0] for idx in xrange(len(elec_x))]
plt.subplot(247, title='Extracellular\npotentials', xlabel='Time [ms]', ylabel='Normalized', xlim=[9, 18])
[plt.plot(cell.tvec, norm_LFP[idx] / np.max(np.abs(norm_LFP[idx])), c=elec_idx_colors[idx], lw=2) for idx in xrange(len(elec_x))]

LFP_amplitude_decay = np.max(electrode.LFP, axis=1) - np.min(electrode.LFP, axis=1)

plt.subplot(244, title='Extracellular\npotential decay', xlabel='x [$\mu m$]', ylabel='mV')
plt.plot(-elec_x, LFP_amplitude_decay)
plt.grid(True)

plt.subplot(248, title='Extracellular\npotential decay', xlabel='x [$\mu m$]', ylabel='Normalized', ylim=[0, 1.1])
plt.plot(-elec_x, LFP_amplitude_decay / np.max(LFP_amplitude_decay))
plt.grid(True)

plt.savefig('exercise_2_spike.png')

