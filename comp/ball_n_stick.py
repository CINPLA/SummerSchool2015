
import pylab as plt
import neuron
nrn = neuron.h

soma = nrn.Section('soma')
soma.L = 15.  # um
soma.diam = 15.
soma.nseg = 1

dend = nrn.Section('dend')
dend.L = 1000
dend.diam = 2
dend.nseg = int(dend.L / 10)

dend.connect(soma, 1, 0)

for sec in nrn.allsec():
    sec.Ra = 100  # Ohm cm
    sec.cm = 1  # uF / cm2
    sec.insert('hh')

stim = nrn.IClamp(soma(0.5))
stim.delay = 20
stim.dur = 200  # ms
stim.amp = -0.2  # nA

t = nrn.Vector()
t.record(nrn._ref_t)

v = nrn.Vector()
v.record(soma(0.5)._ref_v)

nrn.finitialize(-65)
neuron.run(400)

plt.plot(t, v)
plt.savefig('ball_n_stick.png')