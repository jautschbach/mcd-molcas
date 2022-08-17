#!/usr/bin/env python
#import pandas as pd
import argparse
import numpy as np
from os.path import isfile
from sys import exit


###########################
#
# Parsing the number of SO STates
#
###########################
parser = argparse.ArgumentParser(usage='file.py NStates DEBUT',
                                 description='specify the number of states')
parser.add_argument('NStates', type=int, help='number of Spin-Orbit States')
parser.add_argument('DEBUT', type=int, help='starting')
args = parser.parse_args()
#Check if we have all the arguments
DEBUT=args.DEBUT
NStates=args.NStates
if args.NStates == [0]:
    print('Expecting the number of SpinOrbit States')
    exit()
else:
    print('######################### \n# \n#  Now lets check the data:\n# ED for electric dipole, MD for magnetic dipole and Q for electric quadrupole contributions \n#########################') 

def en_file_to_numpy(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    x = []
    for i in lines[1:]:
        x.append(float(i))
    return np.array(x)

def propangmom_file_to_numpy(nstates, filename):
    mat = np.zeros([nstates,nstates], dtype=np.complex_)
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i in range(1, nstates*nstates+1):
        np.set_printoptions(precision=17)
        sp = lines[i].split()
        mat[int(sp[0])-1, int(sp[1])-1] = complex(-float(sp[3]), float(sp[2]))
    return mat

def prop_file_to_numpy(nstates, filename):
    mat = np.zeros([nstates,nstates], dtype=np.complex_)
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i in range(1, nstates*nstates+1):
        np.set_printoptions(precision=17)
        sp = lines[i].split()
        mat[int(sp[0])-1, int(sp[1])-1] = complex(float(sp[2]), float(sp[3]))
    return mat

if not isfile('energies.txt'):
    exit("I don't see energies.txt file? Run molcas-get-energies.sh <rassi.out>")

# energies
en = en_file_to_numpy('energies.txt')
nstates = len(en)

# dipole moments
dipx = prop_file_to_numpy(nstates, 'dipole-1.txt')
dipy = prop_file_to_numpy(nstates, 'dipole-2.txt')
dipz = prop_file_to_numpy(nstates, 'dipole-3.txt')

# angular momenta
angx = propangmom_file_to_numpy(nstates, 'angmom-1.txt')
angy = propangmom_file_to_numpy(nstates, 'angmom-2.txt')
angz = propangmom_file_to_numpy(nstates, 'angmom-3.txt')

# quadrupole moments
quadxx = prop_file_to_numpy(nstates, 'quadrupole-1.txt')
quadxy = prop_file_to_numpy(nstates, 'quadrupole-2.txt')
quadxz = prop_file_to_numpy(nstates, 'quadrupole-3.txt')
quadyy = prop_file_to_numpy(nstates, 'quadrupole-4.txt')
quadyz = prop_file_to_numpy(nstates, 'quadrupole-5.txt')
quadzz = prop_file_to_numpy(nstates, 'quadrupole-6.txt')

# spins
spinx = prop_file_to_numpy(nstates, 'spin-1.txt')
spiny = propangmom_file_to_numpy(nstates, 'spin-2.txt')
spinz = prop_file_to_numpy(nstates, 'spin-3.txt')

# Some Constants
au2ev = (2.7211386021e1)
au2cm = (2.1947463e5)
ev2au = (1.0 / au2ev)
ge = 2.00231930436182
d_au2cgs = (6.460475024e-36)
#alpha is used to convert Atomic units to cgs
#it is mentioned in Pedersen Chem Phys Lett (1995) 246 1
alpha = (4.7144364e2)
c_au = (137.03599914)

print('#   E(eV)   |    E(cm)    |   D_ED(CGS)   | D_ED_MD(CGS) |    fED     |    fMD     |    fEQ     |    R(CGS)    |    g_ED    |  g_ED_MD   |  g_ED_MD_Q')
print(160 * '-')

final_state = NStates
debut_state = DEBUT
for i in range(debut_state+1,final_state):
    # energies
    E_au = en[i] - en[debut_state]
    E_ev = E_au * au2ev
    E_cm = E_au * au2cm

    #################################################################
    # Electric Dipole and its oscillator strength
    #
    # Oscillator strength calculated according to Eq 19 of
    # S. DeBeer et al. Inorganica Chimica Acta 361 (2008) 965
    #
    #################################################################
    Di = np.array([dipx[debut_state][i], dipy[debut_state][i], dipz[debut_state][i]])
    f = (2.0/3.0) * E_au * ( np.sum( Di * Di.conjugate() ) ).real
    D_au = (3.0 * f) / (2.0 * E_au)
    D_cgs = D_au * d_au2cgs

    #################################################################
    # Magnetic Dipole and its oscillator strength
    #
    # Oscillator strength calculated according to Eq 20 of
    # S. DeBeer et al. Inorganica Chimica Acta 361 (2008) 965
    #
    #  The operator is 0.5*(L + ge*S)
    #
    #################################################################
    # --> get amgmom component
    Mi = np.array([angx[debut_state][i], angy[debut_state][i], angz[debut_state][i]])
    # --> include the spin component
    Mi += ge * np.array([spinx[debut_state][i], spiny[debut_state][i], spinz[debut_state][i]])
    # --> include the prefactor of 1/2
    Mi *= 0.5
    # --> get the oscilator strenghts
    fMD = ((2 * E_au) / (3)) * ((1/c_au)**2) * ( np.sum(Mi * Mi.conjugate()) ).real
    MD_au = (3.0 * fMD * (c_au**2)) / (2.0 * E_au)
    MD_cgs = MD_au * d_au2cgs

    #################################################################
    # Electric Quadrupole and its oscillator strength
    #
    # Oscillator strength calculated according to Eq 21 of
    # S. DeBeer et al. Inorganica Chimica Acta 361 (2008) 965
    #
    #################################################################
    ## primitive quadrupole tensor
    qab = np.array([
            [quadxx[debut_state][i], quadxy[debut_state][i], quadxz[debut_state][i]],
            [quadxy[debut_state][i], quadyy[debut_state][i], quadyz[debut_state][i]],
            [quadxz[debut_state][i], quadyz[debut_state][i], quadzz[debut_state][i]],
    ])
    # form traceless quadrupole tensor
    trace = np.trace(qab)
    traceless_qab = np.zeros([3, 3], dtype=np.complex_)
    for alpha_loop in range(3):
        for beta_loop in range(3):
            if alpha_loop == beta_loop:
                traceless_qab[alpha_loop][beta_loop] = (3.0 * qab[alpha_loop][beta_loop] - trace) / 2.0
            else:
                traceless_qab[alpha_loop][beta_loop] = (3.0 * qab[alpha_loop][beta_loop]) / 2.0
    
    # compute f
    trace = np.square(np.trace(traceless_qab)) / 3.0
    val = np.sum(traceless_qab * traceless_qab.conjugate() - trace)
    fEQ = (E_au ** 3) / (20 * c_au ** 2) * val.real
    
    #################################################################
    # Rotatory strength for the isotropic case
    #
    # Rotatory strength calculated according to Eq 24 of
    # Richardson J Chem Phys (1982) 76, 1595
    #
    #################################################################
    R_au = ( np.sum( Mi * Di.conjugate() ) ).imag
    R_cgs = R_au * alpha *10 ** (-40)
    
    #################################################################
    #
    #  Transition Dipole Elements
    #
    #################################################################
    # Containing electric and magnetic dipole contributions
    mag_approx_f = f + fMD
    mag_approx_D_au = (3.0 * mag_approx_f) / (2.0 * E_au)
    mag_approx_D_cgs = mag_approx_D_au * d_au2cgs
    # Containing electric, magnetic and electric quadrupole contributions 
    quad_approx_f = f + fMD + fEQ
    quad_approx_D_au = (3.0 * quad_approx_f) / (2.0 * E_au)
    quad_approx_D_cgs = quad_approx_D_au * d_au2cgs

    #################################################################
    #
    # Disymmetry Factor
    #
    #################################################################
    g_ED = 4 * R_cgs / D_cgs
    g_ED_MD = 4 * R_cgs / mag_approx_D_cgs 
    g_ED_MD_Q = 4 * R_cgs / quad_approx_D_cgs

    #################################################################
    #
    # Lest's print all the results now
    #
    #################################################################
    print('{:11.5f} | {:11.5f} | {:12.5e}  | {:12.5e} | {:1.4e} | {:1.4e} | {:1.4e} | {:12.5e} | {:10.5f} | {:10.5f} | {:10.5f}'.format(E_ev, E_cm, D_cgs, mag_approx_D_cgs, f, fMD, fEQ, R_cgs, g_ED, g_ED_MD, g_ED_MD_Q))
    #print(80 * '-')
