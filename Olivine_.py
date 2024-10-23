#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 15:49:56 2021

@author: jonw
"""

import copy

def Olivine_(ol_percent, initial, elements):
    """
    

    Parameters
    ----------
    ol_percent : integer
        percent of Olivine in liquid
    initial : dictionary
        the starting liquid composition - oxide wt%
    elements : list
        major elements of silicate liquid

    Returns
    -------
    silicate_ox : TYPE
        DESCRIPTION.
    ol_D : dictionary
        olivine/liquid partition coefficient
    ol_ats : dictionary
        olivine composition (oxide %)

    """
    # loop over the number of olivine fractionation steps to recalc liquid comps
    silicate_ox = copy.deepcopy(initial)
    
    
    _els = [ 'Si', 'Ti', 'Al', 'Fe', 'Ca',  'Na', 'K', 'Mg']
    ox = [60.0843, 79.8658, 101.9613, 71.8444, 56.0774, 61.9789, 94.1960, 40.3044]
    oxide_mass = {_els: row for _els, row in zip(_els, ox)}
    


    for ni in range(ol_percent):

        
           #calculate the Mg# of liquid and olivine:

        # mol props

        XMg = silicate_ox['Mg']/oxide_mass['Mg']
        XFe = silicate_ox['Fe']/oxide_mass['Fe']
        

        kd = 1/0.32
        Mg_no_ol = kd*(XMg/XFe)/(1.+kd*(XMg/XFe))
        

        # calc composition of Ol in equilibrium with the liquid

        ol_ats = []

        ol_si_moles = 60.0843
        ol_mg_moles = Mg_no_ol*2.*oxide_mass['Mg']

        ol_fe_moles = (1-Mg_no_ol)*oxide_mass['Fe']

    
        total_olivine_moles = ol_si_moles+ol_mg_moles+ol_fe_moles
        
        ol_mg = 100*ol_mg_moles / total_olivine_moles
        ol_si = 100*ol_si_moles / total_olivine_moles
        ol_fe = 100*ol_fe_moles / total_olivine_moles
        ol_wts = {'Si':ol_si, 'Mg': ol_mg, 'Fe': ol_fe}
        

        #  modify olivine elements in melt
        for e in range (len(_els)):

            g = _els[e]
            

            if (g == 'Si' or g =='Fe' or g =='Mg'):

                # only change the values for these elements

                silicate_ox[g] = silicate_ox[g]-ol_wts[g]*0.01

                continue

        # normalise all liquid els back to 100 - pain in the arse works
        factor = 100./sum(silicate_ox.values())
     

        
        for h in range (len(_els)):

            j = _els[h]

            silicate_ox[j]= factor * silicate_ox[j]
        
   
    
    
    ol_D = dict.fromkeys(elements, 1e-3)
    for key in ol_D:
        if key == 'Ni':
            ol_D[key] = 3.346* (ol_mg)/(silicate_ox['Mg'])- 3.665
        elif key == 'Co':
            ol_D[key] =  0.786* (ol_mg)/(silicate_ox['Mg']) - .385
        else:
            continue
        
    


    return  silicate_ox, ol_D, ol_ats


def Silicate_moles(sil_wts, oxsum, silicate_els, at_el_to_oxide, atdata, valence):
    
    """
    silicate - silicate wt% (sil_wts)
    oxsum - sum of oxides (sum of sil_ox)
    silicate_els - list of silicate elements, not including Oxygen
    """
    
    oxides = dict.fromkeys(silicate_els, 0.)
    el_moles = dict.fromkeys(silicate_els, 0.)  
    el_prop = dict.fromkeys(silicate_els, 0.)
    
    ox_prop = dict.fromkeys(silicate_els, 0.)
    ox_moles = dict.fromkeys(silicate_els, 0.)

    el_sum = 0.
    el_total_prop = 0.
    
    ox_total_prop = 0.
    ox_sum = 0.
    

    for i in range(len(silicate_els)):
        # element mole fracs
        e = silicate_els[i]
        
        el_sum = el_sum + sil_wts[e]  # count up elements, for Oxygen by dif
        
        el_prop[e] = sil_wts[e]/atdata[e]
        el_total_prop = el_total_prop + el_prop[e]  # mol prop for wts
        
        oxides[e] = sil_wts[e]/at_el_to_oxide[e]
        
        if valence[e] == 1:
            ox_prop[e] = oxides[e]/(2*atdata[e]+atdata['O'])
            
        elif valence[e] == 2:
            ox_prop[e] = oxides[e]/(atdata[e]+atdata['O'])
        elif valence[e] == 3:
            ox_prop[e] = oxides[e]/(2*atdata[e]+3*atdata['O'])
        elif valence[e] == 4:
             ox_prop[e] = oxides[e]/(atdata[e]+2*atdata['O'])
            
        ox_total_prop = ox_total_prop + ox_prop[e]
        
      
    

    oxygen = oxsum - el_sum
    
    el_prop['O'] = oxygen/atdata['O']
    el_total_prop = el_total_prop + el_prop['O']
    
    
    for j in range(len(silicate_els)):
        f = silicate_els[j]
        el_moles[f] = el_prop[f]/el_total_prop  # element mole fractions
        ox_moles[f] = ox_prop[f]/ox_total_prop

    return el_moles, ox_moles


