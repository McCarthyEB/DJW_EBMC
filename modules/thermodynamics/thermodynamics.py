# -*- coding: utf-8 -*-
"""
Routines for carrying out thermodynamic calculations

Started on Thu Mar 25  2021

@author: dave
"""
import csv
import sys
import os

from physical_constants import *

# Will use numpy for arrays

import numpy as np
from scipy import special

from ase.io import read, write, Trajectory

from ase import Atom
from ase import Atoms
#
# Work out the partition function and average adsorption energy for a particular 
# list of energies, Eads[], and degeneracies, degen[], 
# at a particular temperature, temp.
# Function also returns the contrib array giving the contribution of each member of
# the energy list to the partition function. This array is normalised.
# report is a logical flag to ask for printing.
#
def calc_Eave(degen, Eads, temp, report):
#
    Emin=np.amin(Eads)
    if report:
       print("Emin= %10.6f" % Emin)
#
# Pass an array containing the contribution of each configuration to the average.
#
    contrib=[]
#
    kt_eV= k_boltz*temp/e_charge
    if report:
       print("kt_eV = %10.6f" % kt_eV)
    q=0
    Esum=0
    for iii in range(0,len(Eads)):        
# 
        Erel=Eads[iii]-Emin
#
        boltz=np.exp(-Erel/kt_eV)
        contrib.append(degen[iii]*boltz)
#
        contrib_E=contrib[iii]*Erel
        q= q + contrib[iii]
        if report:
           print("Eads = %10.6f, degen = %d, E_rel = %10.6f contrib: %10.6f q sum: %10.6f" \
                                    % (Eads[iii], degen[iii], Erel, contrib[iii], q) )
#
        Esum= Esum + contrib_E
#
    Eexc=Esum/q
    Eav=Emin+Eexc
#
    if report:
         print("Minimum energy: %10.6f eV average: %10.6f eV." \
                   % (Emin, Eav))
#
# Normalise the contributions array
#
    for iii in range(0,len(Eads)):        
       contrib[iii] = contrib[iii]/q
    
    return(Eav, Emin, q, contrib)
# 
# Work out the zpe internal energy and vibrational entropy for a set of frequencies at given temperature
# contrib is the array of contributions for the configurations.
# is_2Dgas allows calculation of modes only in z-direction. The modes parrallel to the surface are
# treated as a gas so no vibrations.
# Assume that Eads is the calculated adsorption energy list and does not have ZPE difference already included
# mol_vib is the gas phase reference vibrational frequency.                                                  
#
def freq_list_ZPE_U_S(Eads, freq_set, nads, contrib, barrier, temp, mol_vib, use_hindered, is_2Dgas, give_zpe):
# 
    ncons=len(Eads)
    uuu=[]
    uhind=[]
    sss=[]
    zzzpe=[]
    contrib_freq=[]
    Hads=[]
#
# 5kT/2 for translation and rotational degrees of freedom
# additional kT for the pV term in the enthalpy
#
    hgas = 7.0*temp*k_boltz/(2*e_charge)
#
# 100 required as  mol_vib waveno in cm-1.
    vib_ene_ev= mol_vib*h_planck*c_light*100.0/e_charge
#
    zpe_gas = 0.5*vib_ene_ev
#
    if give_zpe:
       print("In freq_list_ZPE_U_S Eads has length: %d, nads: %d, len freq_set: %d" % (len(Eads), nads, len(freq_set)))
       if use_hindered:
          print("Will include hindered modes....")
       else:
          print("Will only treat modes harmonically")
#
       if is_2Dgas:
          print("Will treat xy as 2D gas........")
       else:
          print("This is not a 2D gas")
#
# Assume highest 1/3rd of frequencies are perpendicular to surface for each configuration
# configurations each have 3*nads frequencies
#
    nf_per_con=int(3*nads)
    two_thirds_nf=2*nf_per_con/3
    icon=0
    which_con=0
    zpe_config=0.0
    uuu_config=0.0
    uhind_config=0.0
    Hads.append(Eads[0])
#    print("Calculating vib freq at %4.1f K nf_per_con = %d two thirds: %d" % (temp, nf_per_con, two_thirds_nf))
#
# nfreqs is the total number of frequencies for all configurations
#
    nfreqs=len(freq_set)

    if give_zpe:
        print("Have %d frequencies altogether...." % nfreqs)
#
# Now that freq_set contains all the vibrations for the whole set of configurations
# only need to look over freq_set itself to get the vibrational contributions to U and S.
#
    for ifreq in range(0,nfreqs):
        is_xy = (icon < two_thirds_nf)
        if (use_hindered or is_2Dgas):
           is_hindered = is_xy
        else:
           is_hindered = False
#
# U_and_S_vibs_hindtrans returns the contribution to each term for the single vibration indexed by ifreq
# note that for hindered translations the uhind is added into ucontrib by this routine and so its values is
# only returned for comparison purposes.
#
        zpe_contrib,ucontrib,uhindcontrib,scontrib=U_and_S_vibs_hindtrans(freq_set[ifreq], temp, \
                                                                          barrier, is_hindered, is_2Dgas)
#
#  contrib_freq provides the scale factor for individual frequencies so that they are scaled to
#  a per atom basis.
#
        if not is_2Dgas:
            contrib_freq.append(contrib[which_con]/float(nads))
        else:
            if is_xy:
              contrib_freq.append(0) 
            else:
              contrib_freq.append(contrib[which_con]/float(nads)) 
#
# contributions of vibration ifreq to entropy and internal energy with weighting for configuration
#
        sss.append(contrib_freq[ifreq]*scontrib)
        uuu.append(contrib_freq[ifreq]*ucontrib)
        uhind.append(contrib_freq[ifreq]*uhindcontrib)
        zzzpe.append(contrib_freq[ifreq]*zpe_contrib)
#
# accumulate the ZPE contributions for the configurations in eV per H2 molecule units
# The contrib factor will be taken into account when Hav is formed so here just need to
# ensure the scale factor is correct for a per H2 molecule unit.
# If this is_2Dgas only the mode perpendicular to the surface is to be included.
# for lattice gas and hindered all three modes are used so contribute.
#
        if not is_2Dgas:
            freq_fac= 2.0/float(nads)
        else:
            if is_xy:
              freq_fac = 0.0
            else:
              freq_fac= 2.0/float(nads)
#
        zpe_config=zpe_config + freq_fac*zpe_contrib
        uuu_config=uuu_config + freq_fac*ucontrib
        uhind_config=uhind_config + freq_fac*uhindcontrib
#
#        if is_hindered:
#           print("freq %d %10.6f cm-1 icon %d is hindered, ZPE: %10.6f eV, U contrib %10.6f eV, temp*scontrib: %10.6f eV weight %10.6f" % \
#                                    (ifreq, freq_set[ifreq], icon, zpe_contrib, ucontrib, temp*scontrib/e_charge, contrib_freq[ifreq]))
#        else:
#           print("freq %d %10.6f cm-1 icon %d is free, ZPE: %10.6f eV,  U contrib %10.6f eV, temp*scontrib: %10.6f eV weight %10.6f" % \
#                                    (ifreq, freq_set[ifreq], icon, zpe_contrib, ucontrib, temp*scontrib/e_charge, contrib_freq[ifreq]))
        icon=icon+1
#
# Done frequencies for this configuration so can update Hads
#
        if icon == nf_per_con:
#
           if (give_zpe):
              print("icon = %d " % icon)
              print("Config %d has Eads %10.6f eV per H2" % (which_con, Eads[which_con]))
              print("Config %d has Hads first term %10.6f eV per H2" % (which_con, Hads[which_con]))
              print("Config %d ZPE per H2 : %10.6f eV zpe_gas ref: %10.6f eV" % (which_con, zpe_config, zpe_gas))
              print("Config %d uuu_config: %10.6f eV uhind_config: %10.6f eV" % (which_con, uuu_config, uhind_config))
              print("Will use references zpe_gas: %10.6f hgas: %10.6f\n" % ( zpe_gas, hgas))
#
           Hads[which_con]=Hads[which_con] + zpe_config + uuu_config 
           zpe_config=0.0
           uuu_config=0.0
           uhind_config=0.0
           icon=0
           which_con=which_con+1
#           print("which_con: %d ncons: %d len Hads: %d" % ( which_con, ncons, len(Hads)))
           if which_con < ncons:
              Hads.append(Eads[which_con])
#
# These are sums over the full set of frequencies but the contrib arrays ensure they
# are scaled per atom ( 3 modes ).
#
    Svibs=np.sum(sss)/e_charge
    Uvibs=np.sum(uuu)
    Uhinds=np.sum(uhind)
    ZPEvibs=np.sum(zzzpe)
#
# Take the hgas and zpe_gas term off Hads
#
    for iii in range(0,len(Hads)):
       Hads[iii]=Hads[iii] - zpe_gas - hgas
#
    return(ZPEvibs, Uvibs, Uhinds, Svibs, contrib_freq, Hads)
#
# rotational partition function for H2 using summation approach
# This breaks the sum into even and odd terms to account for para and ortho H2.    
#
def q_rot_summed_H2(rot_con, temp):
#
# acc defines the smallest contribution that should be included, 
# i.e. the cutoff for the summation.
#
    acc=1E-6
    B_ev=rot_con*100.0*h_planck*c_light/e_charge
    kt_ev=k_boltz*temp/e_charge
#
# Even terms
#
    j=2
    q_sum_even=0.0
    contrib=1.0
    while contrib > acc:
       q_sum_even= q_sum_even+contrib
       jterm=j*(j+1)
       contrib= (2*j+1)*np.exp(-jterm*B_ev/kt_ev)
       j=j+2
#
# Odd terms
#
    j=-1
    q_sum_odd=0.0
    contrib=1.0
    while contrib > acc:
       j=j+2
       jterm=j*(j+1)
       contrib= (2*j+1)*np.exp(-jterm*B_ev/kt_ev)
       q_sum_odd= q_sum_odd+contrib
#
    q_rot= q_sum_even + 3*q_sum_odd 
#
    return q_rot
#
# Work out gas phase partition function and absolute entropy at this T and p
# Written for diatomics with mass in relative atomic mass units
# vib_freq and rotational constant in cm-1 units
# symfac is the rotational symmetry factor, 2 for O2 etc 1 for HCl etc
#        for the case of H2 use the q_rot_summed function, if using the
#        analytic form set symfac = 0.5 to account for average nuclear spin degeneracy
# temperature should be in Kelvin and pressure in Pa.
# spin is the spin state of the molecule, 1 for H2 singlet, 3 for O2 triplet
#
def Gas_q_and_S(mass, vib_waveno, rot_con, symfac, spin, temp, press):
#
# kt in SI units and in eV units 
    kt_si=k_boltz*temp
    kt_ev = kt_si/e_charge
#
# translational contribution
#
    c=2.0*np.pi*amu*mass*k_boltz/h_planck
    c=c/h_planck
    q_trans_prefactor = R_gas*np.power(c,1.5)    
#
# translational partition function
#
    q_trans= q_trans_prefactor*np.power(temp,2.5)/press
#
# get translational contribution in eV units
# Note that np.log() is the natural logarithm
#
# build entropy for mole of H2 from 3R/2 + R ( ln qt - ln NA + 1) as the translational contrib
#
    trans_contrib= R_gas*(3/2 + np.log(q_trans) - np.log(n_avo) +1 )
#
# Rotational contribution, 100 required as rot_con in cm-1.
#
#    B_si=rot_con*100.0*h_planck*c_light
#    q_rot= kt_si/(symfac*B_si)
# Use summation version with ortho/para hydrogen contributions
    q_rot=q_rot_summed_H2(rot_con, temp)
#
    rot_contrib= R_gas*(1+np.log(q_rot))    #### 1 here for U-U(0)/T term 
                                            #### assuming the equipartition for rotations
#
# Vibrational contribution, 100 required as vib_waveno in cm-1.
# also calculate gas phase ZPE in eV units
#
    vib_ene= vib_waveno*h_planck*c_light*100.0
    vib_ene_ev = vib_ene/e_charge
    zpe=0.5*vib_ene_ev
    earg= -vib_ene_ev/kt_ev
    q_vib= 1/(1-np.exp(earg))
#
    vib_contrib= R_gas*(-earg/(np.exp(-earg)-1) + np.log(q_vib))
#
# Spin contribution
#
    spin_contrib = R_gas*np.log(spin)
#
#    if (temp > 40  and temp < 350):
#       beta=1.0/(k_boltz*temp)
#       lambdaw=np.sqrt(1.0/(temp*c))
#       print("Data in Gas_q_and_S with temp= %10.6f K pressure = %10.6e Pa" % ( temp, press ))
#       print("beta: %10.6e J-1, lambda: %10.6e m, c: %10.6e" % (beta, lambdaw, c))
#       print("Partition functions: q_trans = %10.6e, q_rot = %10.6f, q_vib = %10.6f" \
#              % ( q_trans, q_rot, q_vib ))
#       print("Absolute Entropy for gas contribs (eV/K) trans: %10.6f, rot: %10.6f, vib: %10.6f, spin: %10.6f" \
#              % ( trans_contrib/(n_avo*e_charge), rot_contrib/(n_avo*e_charge), \
#                  vib_contrib/(n_avo*e_charge), spin_contrib//(n_avo*e_charge)))
#
    S_gas = (trans_contrib + rot_contrib + vib_contrib + spin_contrib)
# convert S_gas to eV K-1
    S_gas = S_gas/(n_avo*e_charge)
# also send back translational contributuion
    S_trans= trans_contrib/(n_avo*e_charge)
#
    q_all = q_trans*q_rot*q_vib*spin
    
    return(q_all, S_gas, S_trans, zpe)
#
# Work out chemical potential shift for this T and p
# Written for diatomics with mass in relative atomic mass units
# vib_freq and rotational constant in cm-1 units
# symfac is the rotational symmetry factor, 2 for H2, O2 etc 1 for HCl etc
# temperature should be in Kelvin and pressure in Pa.
# spin is the spin state of the molecule, 1 for H2 singlet, 3 for O2 triplet
#
def del_chem_pot(mass, vib_waveno, rot_con, symfac, spin, temp, press):
#
# translational contribution
#
    c=2.0*np.pi*amu*mass/h_planck
    c=c/h_planck
    q_trans_prefactor = np.power(c,1.5)    
#
# get kT in SI and eV units
#
    kt= k_boltz*temp
    kt_ev= k_boltz*temp/e_charge
    argt= q_trans_prefactor*np.power(kt,2.5)/press
#
# get translational contribution in eV units
# Note that np.log() is the natural logarithm
# putting factor of 0.5 in for each contrib would mean 
# we are calculating the molecular 
# partition function but per atom.
#
    trans_contrib= -kt_ev*np.log(argt)
#
# Rotational contribution, 100 required as rot_con in cm-1.
#  
#    B_si=rot_con*100.0*h_planck*c_light
#    arg= kt/(symfac*B_si)
# Use summation version with ortho/para hydrogen contributions
    q_rot=q_rot_summed_H2(rot_con, temp)
    rot_contrib= -kt_ev*np.log(q_rot)
#
# Vibrational contribution to gas phase chemical potential, 
# 100 required as vib_waveno in cm-1.
#
    vib_ene= vib_waveno*h_planck*c_light*100.0
    earg= -vib_ene/kt
    arg= 1/(1-np.exp(earg))
    vib_contrib= -kt_ev*np.log(arg)
#
# Spin contribution
#
    spin_contrib = -kt_ev*np.log(spin)
#
    del_mu = (trans_contrib + rot_contrib + vib_contrib + spin_contrib)
#    if (temp > 250 and temp < 350):
#      mass_term=2.0*np.pi*amu*mass
#       print("In del_chem_pot: mass_term = %10.6e, c = %10.6e,  q_trans_prefactor = %10.6e, argt = %10.6e" % \
#                                     (mass_term, c, q_trans_prefactor, argt))
#       print("kt_ev = %10.6e" % kt_ev)
#       print("del_mu contribs at %5.2f K trans: %10.6e eV, rot: %10.6e eV, vib: %10.6e eV, spin: %10.6e eV" \
#              % ( temp, trans_contrib, rot_contrib, vib_contrib, spin_contrib))
    
    return(del_mu)
#
# Calculate internal energy contribution and entropy of adsorbate vibrational states 
# using Campbell's hindered translator
# equations from J. Phys. Chem. C, 120, 9719-9731, (2016).
#
# barrier is the barrier to diffusion in eV units
#
def U_and_S_vibs_hindtrans(vib_waveno, temp, barrier, is_hindered, is_2Dgas):
#
# Catch zero or negative vibs and return zero
   Uhind = 0.0
   if vib_waveno < 10:
      return (0,0,0,0)
#
# 100 required as vib_waveno in cm-1.
   vib_ene_si= vib_waveno*h_planck*c_light*100.0
   vib_ene_ev= vib_ene_si/e_charge
#
   zpe = 0.5*vib_ene_ev
#
# kT in SI and in eV units
#
   kt_si = k_boltz*temp
   kt_ev = kt_si/e_charge
#
# dimensionless temperature and its inverse
#
   temp_nodim = kt_si/vib_ene_si
   temp_nodim_inv = 1.0/temp_nodim
#
# terms for all modes
# SHO contribution to internal energy and entropy
# is_hindered can be taken as indicating modes parallel to surface, for 2Dgas these are 
# found separately based on the 3D gas translational result.
#
   if (is_hindered and is_2Dgas):
      Uho = 0.0
      Sho = 0.0
   else:
      Uho = kt_ev*temp_nodim_inv/(np.exp(temp_nodim_inv)-1)
#
      Sho = k_boltz*temp_nodim_inv/(np.exp(temp_nodim_inv)-1 ) 
      Sho = Sho - k_boltz*( np.log(1.0 - np.exp( -temp_nodim_inv )) )
#
# Additional terms for hindered translations
#
   if (is_hindered and not is_2Dgas):
#
      ratio = barrier/vib_ene_ev
#      print("DEBUG: Hindered mode, barrier %10.6f eV, vib_ene: %10.6f eV, ratio %10.6f" % (barrier, vib_ene_ev, ratio))
#
      arg= 0.5*ratio*temp_nodim_inv
      fac= np.sqrt(ratio*temp_nodim_inv*np.pi)
#
      i0 = special.i0(arg)
      i1 = special.i1(arg)
#
      Uhind = -1.0*kt_ev*(0.5 + temp_nodim_inv/(2+16*ratio)-arg*(1-i1/i0))
      Uho= Uho + Uhind
#
      Sho= Sho - k_boltz*(0.5 + arg*i1/i0 - np.log(fac*i0) )
#
# Note Uho is returned in eV units, Sho is returned in SI units J K-1
# Note for hindered translations Uho contains the hindered contribution so
#      Uhind only returned for comparison.
#
   return (zpe,Uho,Uhind,Sho)
