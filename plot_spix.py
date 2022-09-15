import numpy as np
import numpy.ma as ma


def plot_spix(img1,img2,freq1,freq2,noise1,noise2,sigma=3):
    spix1 = img1*(img1 > noise1*sigma) #replaces indices where condition is not met with 0
    spix2 = img2*(img2 > noise1*sigma)
    a = np.log10(spix1/spix2)/np.log10(freq1/freq2)

    level10 = noise1*sigma
    lev1=[]
    level20 = noise2*sigma
    lev2=[]

    for i in range(0,10):
        lev1.append(level10*2**i)
        lev2.append(level20*2**i)

    f,ax = plt.subplots()

    cset = ax.contour(img1,lev1,inline=1,colors=['grey'], extent=extent, aspec=1.0)
    im = ax.imshow(a,origin='bottom',extend= extent, vmin=vmin, vmax=vmax)

    ax.axis('scaled')
  ax.set_xlabel('Right Ascension [mas]',fontsize=15)
  ax.set_ylabel('Relative Declination [mas]',fontsize=15)
    ax.set_xlim()
    ax.set_ylim()
    ax.minorticks_on()
    ax.tick_params('both', length=8, width=2, which='major')
    f.tick_params(axis='both',which='both',direction='in', labelsize=13)

    cb = plt.colorbar(im,cax=cax,cmap='jet')
    cb.ax.tick_params(labelsize=13, width=2)
    cb.set_label(r'$\alpha$',fontsize=15)

    plt.savefig('spectral_index_between_{.1f}_{.1f}.pdf'.format(freq1,freq2),bbox_inches='tight')
    plt.close('all')
 
