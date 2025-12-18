import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.interpolate
import os


def plot_powerscan(date='2025-12-18', run=3, lineout_skip=8, ignore_rows=0, xmin=1000, xmax=2800, ymin=-70, ymax=-10, cmin_offset=5, plot_powers=8):

    fig, axs = plt.subplots(2, 1, figsize=(8, 5), tight_layout=True, sharex=True)
    ax0, ax1 = axs

    filename = 'data/' + date + '/run_%04i.txt' % run
    
    with open(filename, 'r') as f:
        l = f.readline()    
        osa_type = l.split()[-1]
    
    skiprows = 1
    with open(filename, 'r') as f:
        searching = True
        while searching:
            l = f.readline()
            if "Data:" in l:
                searching=False
            elif "Powers:" in l:
                pows = np.array([float(x) for x in l.split()[1:]])
            else:
                skiprows += 1
    
        wls_str = f.readline()
        wls = np.array([float(x) for x in wls_str.split()])
        
        spectra = np.ones((pows.size, wls.size)) * np.nan
        lines = f.readlines()
        count = 0
        processing = True
        print(len(lines))
        for l in lines:
            print('%i / %i  - power: %.1f mW'%(count+1, pows.size, pows[count]))
            splits = np.array([float(x) for x in l.split()])
            if splits.size == wls.size + 3:
                spectra[count] = splits[3:]
            else:
                break
            count+=1
        # skiprows += 1
        
    # faster loading method that isn't robust to incomplete files
    # data = np.loadtxt(filename, skiprows=skiprows)
    # nums = data[:, 0]
    # pows = data[:, 1]
    # # degs = data[:, 2]
    # spectra = data[:, 3:]

    if xmin is None:
        xmin = np.nanmin(wls)
    if xmax is None:
        xmax = np.nanmax(wls)
    if ymin is None:
        ymin = np.nanmin(spectra) - 10
    if ymax is None:
        ymax = np.nanmax(spectra) + 10
        
    if len(pows) < ignore_rows-2:
        raise ValueError('ignore_rows cannot be larger than the length of the data minus 2.\nIgnore rows: %i versus data length: %i'%(ignore_rows, len(data)))
    ax1.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.95, 0.05, int((len(pows)-ignore_rows)/lineout_skip + 1) )))
    extent = (wls.min(), wls.max(), pows.min(), pows.max())
    
    for n in np.arange(1, pows.size-ignore_rows, lineout_skip):
        ax1.plot(wls, spectra[n], label='%.2f mW'%pows[n])
    
    # ax1.plot(wls, spectra[-1], label='%.2f mW'%pows[-1])

    print(spectra)

    if osa_type == 'thorlabs':
        # the Thorlabs OSA records with even frequency spacing.
        # Re-interpolate to even wavelength spacing
        new_wls = np.linspace(xmin, xmax, 1000)
        extent = (xmin, xmax, pows.min(), pows.max())

        new_spectra = np.zeros((pows.size, new_wls.size))
        for n, spec in enumerate(spectra):
            interp = scipy.interpolate.interp1d(wls, spec)
            new_spectra[n] = interp(new_wls)
        
        spectra = np.array(new_spectra)
        wls = new_wls
            

        
        
    axs[0].imshow(spectra, extent=extent, aspect='auto', cmap='jet', interpolation='nearest', clim=(ymin + cmin_offset, ymax))

    axs[0].set_ylabel('Power (mW)')

    axs[0].axvline(1263, color='w', ls='dashed', lw=1)

    axs[1].grid(alpha=0.3)
    axs[1].set_xlabel('Wavelength (nm)')
    axs[1].set_ylabel('Spectral flux (dBm/nm)')
    #ax1.legend(labelcolor='linecolor', fontsize=8, bbox_to_anchor = (1.1,1))
    axs[1].legend(labelcolor='linecolor', fontsize=8)
    axs[1].set_ylim(ymin, ymax)
    axs[1].set_xlim(xmin, xmax)
    
    # for ax in axs:
        # ax.axvline(1762, color='k', lw=1, ls='dashed')

    axs[0].set_title(date + ' - Run %04i' % run )
    if not os.path.exists('plots'):
        os.mkdir('plots')
    fig.savefig('plots/' + date + ' - Run %04i' % run + '.png', dpi=200)
    
    if plot_powers is not None:

        fig1, ax2 = plt.subplots(1, 1, figsize=(
            6, 2.5), tight_layout=True)
    
        def plot_power(idx):
            l, = ax2.plot(wls, spectra[idx], label='Power: %.2f mW' % (pows[idx]))
            # color = plt.getp(l, 'color')
            #ax1.plot(wls, spectra[idx], color=color)
        
        def isiterable(p_object):
            try:
                it = iter(p_object)
            except TypeError:
                return False
            return True

    
        if isiterable(plot_powers):
            for pp in plot_powers:
                idx = np.argmin(np.abs(pows - pp))
                plot_power(idx)
        else:
            idx = np.argmin(np.abs(pows - plot_powers))
            plot_power(idx)
        
        ax2.set_xticks(np.arange(xmin, xmax+200, 200))
        ax2.set_xticks(np.arange(xmin, xmax+100, 100), minor=True)
        ax2.legend(labelcolor='linecolor')
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Spectral flux (dBm/nm)')
        #ax2.axvline(780, color = 'black', alpha = 0.25, linestyle = '--')
        ax2.set_ylim([ymin, ymax])
        ax2.set_xlim([xmin, xmax])
        ax2.grid(alpha=0.2, color='k')
        #ax2.legend(labelcolor='linecolor', bbox_to_anchor = (1,1.5))
        fig1.savefig('Spectrum.png', dpi=400, bbox_inches = 'tight')


plot_powerscan()


plt.show()
