import scipy as sp
import pylab as plt
from scipy.integrate import odeint

## Full Hodgkin-Huxley Model

def simulate(params):
    # External current
    t = params['t']
    if params['protocol'] == 'clamp':
        print 'Applying voltage clamp protocol'
        if params['prop'] == 'passive':
            a = 10
            b = 10
        elif params['prop'] == 'active':
            a = 40
            b = 40
        def v_clamp(t): return -80+a*(t>100)-b*(t>150)  #set voltage clamp protocol
        def I_inj(t, V): # apply voltage clamp.
            return  (v_clamp(t)-V)/(params['R_clamp']*params['C_m'])
    elif params['protocol'] == 'pace':
        print 'Apply pacing stimulus'
        def I_inj(t,V): # apply pacing stimulus
            return 10*(t>100) - 10*(t>200) + 35*(t>300) #10 uA/cm^2 100ms-200ms and 35 uA/cm^2 for t>300
    else:
        print 'Protocol not recognized'

    # Membrane currents (in uA/cm^2)
    #  Sodium (Na = element name)
    def I_Na(V,m,h):return params['g_Na'] * m**3 * h * (V - params['E_Na'])
    #  Potassium (K = element name)
    def I_K(V, n):  return params['g_K']  * n**4     * (V - params['E_K'])
    #  Leak
    def I_L(V):     return params['g_L']             * (V - params['E_L'])
    # Channel gating kinetics
    # Functions of membrane voltage
    def alpha_m(V): return 0.1*(V+40.0)/(1.0 - sp.exp(-(V+40.0) / 10.0))
    def beta_m(V):  return 4.0*sp.exp(-(V+65.0) / 18.0)
    def alpha_h(V): return 0.07*sp.exp(-(V+65.0) / 20.0)
    def beta_h(V):  return 1.0/(1.0 + sp.exp(-(V+35.0) / 10.0))
    def alpha_n(V): return 0.01*(V+55.0)/(1.0 - sp.exp(-(V+55.0) / 10.0))
    def beta_n(V):  return 0.125*sp.exp(-(V+65) / 80.0)
    # Integrate!
    def dALLdt(X, t):
        V, m, h, n = X
        #calculate membrane potential & activation variables
        dVdt = (I_inj(t, V) - I_Na(V, m, h) - I_K(V, n) - I_L(V))
        dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
        dhdt = alpha_h(V)*(1.0-h) - beta_h(V)*h
        dndt = alpha_n(V)*(1.0-n) - beta_n(V)*n
        return dVdt, dmdt, dhdt, dndt

    if params['iter_R_clamp']:
        lineColor = sp.linspace(0,0.7,len(params['R_clamps']))
        keytitle = 'R_clamps'
    elif params['iter_g_L']:
        lineColor = sp.linspace(0,0.7,len(params['g_Ls']))
        keytitle = 'g_Ls'
    elif params['iter_g_Na']:
            lineColor = sp.linspace(0,0.7,len(params['g_Nas']))
            keytitle = 'g_Nas'
    else:
        print 'Either "iter_R_clamp" or "iter_g_L" or "iter_g_Na" must be set to 1'
    plt.figure()

    for i,val in enumerate(params[keytitle]):
        if params['iter_R_clamp']:
            params['R_clamp'] = val
        if params['iter_g_L']:
            params['g_L'] = val
        if params['iter_g_Na']:
            params['g_Na'] = val
        X = odeint(dALLdt, [-80, 0.05, 0.6, 0.32], t)
        V = X[:,0]
        m = X[:,1]
        h = X[:,2]
        n = X[:,3]
        ina = I_Na(V,m,h)
        ik = I_K(V, n)
        il = I_L(V)

        if params['prop'] == 'passive':
            if params['close_up']:
                plt.subplot(2,1,1)
                plt.title('Hodgkin-Huxley Neuron')
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.ylim(-82,-68)
                plt.xlim(98,104)
                plt.xticks([])

                plt.subplot(2,1,2)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.xlabel('t (ms)')
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                plt.ylim(-10, 100)
                plt.xlim(98,104)
            elif params['iter_g_L']:
                plt.subplot(3,1,1)
                plt.title('Leak conductance defines passive membrane resistance')
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.ylim(-82,-68)
                plt.xlim(95,170)
                plt.xticks([])

                plt.subplot(3,1,2)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                plt.xlim(95,170)
                plt.xticks([])

                plt.subplot(3,1,3)
                valstr = str(val)
                plt.title('Zoomed')
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.xlabel('t (ms)')
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                plt.ylim(-5, 30)
                plt.xlim(95,170)
            else:
                plt.subplot(2,1,1)
                plt.title('Hodgkin-Huxley Neuron')
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.ylim(-82,-68)
                plt.xlim(95,170)
                plt.xticks([])

                plt.subplot(2,1,2)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.xlabel('t (ms)')
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
         ##       plt.ylim(-10, 50)
                plt.xlim(95,170)
        elif params['prop'] == 'active':
            if params['close_up']:
                plt.subplot(4,1,1)
                plt.title('Hodgkin-Huxley Neuron')
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.xlim(95,160)
                plt.xticks([])

                plt.subplot(4,1,2)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                ##plt.ylim(-60, 60)
                plt.xlim(95,160)
                plt.xticks([])

                plt.subplot(4,1,3)
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.ylim(-42, -38)
                plt.xlim(98,120)
                plt.xticks([])

                plt.subplot(4,1,4)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.xlabel('t (ms)')
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                plt.xlim(98,120)
            else:
                plt.subplot(2,1,1)
                plt.title('Hodgkin-Huxley Neuron')
                plt.plot(t, V, str(lineColor[i]))
                plt.ylabel('V (mV)')
                plt.xlim(95,160)
                plt.xticks([])

                plt.subplot(2,1,2)
                valstr = str(val)
                plt.plot(t, I_inj(t, V), str(lineColor[i]),label = valstr)
                plt.xlabel('t (ms)')
                plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
                if params['iter_R_clamp']:
                    plt.ylim(-2000, 2000)
                plt.xlim(95,160)
        else:
            print 'prop not recognized'
    plt.legend(title = keytitle[0:-1])
