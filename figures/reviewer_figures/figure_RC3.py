
# standard libraries
from math import log10

# nonstandard libraries
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

from madhype.analysis.methods import estimate_match_frequencies 

# library modifications
plt.rcParams["font.family"] = "serif"


def main():

    linewidth = 5
    fs = 18

    w_total = 96
    w_i,w_j = 20,20
    cpw = (100,)

    alpha = 1

    #N = np.array([60,2,2,34]) # observed well counts (w_o, w_i, w_j, w_ij)
    w_ijs = xrange(0,40)

    total_probs,approx_probs = [],[]

    plt.subplots(1,1,figsize = (8,6))

    for w_ij in w_ijs:

        print 'W - ij:',w_ij

        N = np.array([w_total-(w_ij+w_i+w_j),w_i,w_j,w_ij])
        total,_ = get_probs_actual(N,visual=False,display=True)
        
        f = estimate_match_frequencies({
            'w_i':(w_i,),
            'w_j':(w_j,),
            'w_ij':(w_ij,),
            'w_o':(w_total - (w_i + w_j + w_ij),),
            'w_tot':(w_total,),
            'alpha':alpha,
            'cpw':cpw,
            })

        i =  1 - (1 - f['i'])**cpw[0]
        j =  1 - (1 - f['j'])**cpw[0]
        ij = 1 - (1 - f['ij'])**cpw[0]

        print 'Freqs!:',i,j,ij


        approx = ((ij**-1)*multinomial(N, multinomial_probs((i,j,ij))))

        total_probs.append(log10(total))
        approx_probs.append(log10(approx))
        
    plt.plot(w_ijs,total_probs,label = r'$\int\int\int P(w|f)\dotP(f|H)$',linewidth=linewidth)
    plt.plot(w_ijs,approx_probs,label = r'$argmax( P(w|f)\dotP(f|H) )$',linewidth=linewidth)

    leg = plt.legend(fontsize=fs-4)
    leg.get_frame().set_edgecolor('black')

    plt.xlabel(r'$W_{ij}$',fontsize = fs)
    plt.ylabel('$\log(P)$',fontsize = fs)

    plt.tick_params(labelsize = 14)

    plt.savefig('figure_rc3b.png',dpi=300)
    plt.show(block=False)
    raw_input()
    plt.close()


    #N = np.array([66,5,5,20])
    #total,approx = get_probs(N,visual=True,display=False)
    
    
    

def multinomial(N, P):
  n_tot = sum(N)

  coef = 1.0
  n_temp = 0
  for n in N:
    n_temp += n
    coef *= scipy.special.binom(n_temp, n)

  return coef * np.prod(P**N)

def multinomial_probs(freqs):
  fi, fj, fij = freqs
  P = [
      (1-fi)*(1-fj)*(1-fij),
      fi*(1-fj)*(1-fij),
      fj*(1-fi)*(1-fij),
      1 - (1-fi*fj)*(1-fij)
  ]
  return np.array(P)

def get_probs(N,visual=False,display=False):
    count = 32
    lim = (-4,-1./8)

    fs = 18

    f_range = np.logspace(*lim,num=count)
    f_label = [f for f in np.linspace(*lim,num=count)]

    f_index = [i for i,f in enumerate(f_label) if int(f) == float(f)]
    f_label = [f for i,f in enumerate(f_label) if int(f) == float(f)]

    f_i  = f_range
    f_j  = f_range
    f_ij = np.flip(f_range,0)

    probs = np.zeros((count,count))
    max_prob = -999
    max_freqs = (0,0,0)
    xy = (0,0)
    total = 0.

    step = (lim[1] - lim[0]) / (count - 1)

    for a,(i,j) in enumerate(zip(f_i,f_j)):
        for b,ij in enumerate(f_ij):
            probs[a,b] = (ij**-1)*log10(multinomial(N, multinomial_probs((i,j,ij))))
            total += (10**probs[a,b])*(ij*i*(10**(step/2) - 10**(-step/2)))
            if probs[a,b] > max_prob:
                max_prob = probs[a,b]
                max_freqs = ((i,j,ij))
                xy = (b,a)

    if visual == True:
        plt.imshow(probs,vmin=-25, vmax=0)
        plt.colorbar()

        plt.xticks(f_index,np.flip(f_label,0))
        plt.yticks(f_index,f_label)

        plt.ylabel(r'$f_{ij}$',fontsize = fs)
        plt.xlabel(r'$f_{i} = f_{j}$', fontsize = fs)

        plt.scatter(*xy,s = 50,facecolors='none',edgecolors='black')

        plt.savefig('figure_rc3a.png',dpi=300)

        plt.show(block=False)
        raw_input()
        plt.close()

    if display == True:
        print 'Total probability:',total/count
        print 'Max probability:', 10**max_prob
        print 'F-ij:',max_freqs[2]
        print 'F-i:',max_freqs[0]
        print 'F-j:',max_freqs[1]

    return total/count, 10**max_prob

def get_probs_actual(N,visual=False,display=False):
    count = 40
    lim = (-5,-1./8)
    
    lim_ij = (-2.5,-1./16)

    fs = 18

    f_range = np.logspace(*lim,num=count)
    f_ij_range = np.logspace(*lim,num=count)
    f_label = [f for f in np.linspace(*lim,num=count)]

    f_index = [i for i,f in enumerate(f_label) if int(f) == float(f)]
    f_label = [f for i,f in enumerate(f_label) if int(f) == float(f)]

    f_i  = f_range
    f_j  = f_range
    f_ij = np.flip(f_ij_range,0)

    probs = np.zeros((count,count))
    max_prob = -999
    max_freqs = (0,0,0)
    xy = (0,0)
    total = 0.

    step = (lim[1] - lim[0]) / (count - 1)

    for a,i in enumerate(f_i):
        for a,j in enumerate(f_j):
            for b,ij in enumerate(f_ij):

                try:
                    probs = log10((ij**-1)*multinomial(N, multinomial_probs((i,j,ij))))
                    total += (10**probs)*(ij*i*j*(10**(step/2) - 10**(-step/2))**3)
                except ValueError:
                    continue

                if probs > max_prob:

                    max_prob = probs
                    multi_prob = log10(multinomial(N, multinomial_probs((i,j,ij))))
                    max_freqs = ((i,j,ij))

    if visual == True:
        plt.imshow(probs,vmin=-25, vmax=0)
        plt.colorbar()

        plt.xticks(f_index,np.flip(f_label,0))
        plt.yticks(f_index,f_label)

        plt.ylabel(r'$f_{ij}$',fontsize = fs)
        plt.xlabel(r'$f_{i} = f_{j}$', fontsize = fs)

        plt.scatter(*xy,s = 50,facecolors='none',edgecolors='black')

        plt.savefig('figure_rc3a.png',dpi=300)

        plt.show(block=False)
        raw_input()
        plt.close()

    if display == True:
        print 'Total probability:',total/count
        print 'Max probability:', 10**max_prob
        print 'Multinomial probability:', 10**multi_prob
        print 'F-ij:',max_freqs[2]
        print 'F-i:',max_freqs[0]
        print 'F-j:',max_freqs[1]

    return total, 10**multi_prob



if __name__ == "__main__":
    main() 

