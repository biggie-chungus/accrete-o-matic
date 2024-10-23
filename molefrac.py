#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:57:24 2019

@author: jonw
"""
import uncertainties.unumpy as unp



def molefrac_(atdata, element_list,  metalconc, cpres):
    """
    convert wt% of elements in METALS to mole fractions
    Si, Ni, C,  S and O are major elements with others treated as trace.
    This is only relevent if metal comps are calculated using metal/silicate
    partitioning data.  If not, this needs ammending and/or  commented
    out and directly inputting metal compositions
    """

    solutetotal = 0.
    
    total = 0.
    mol_prop_element = dict.fromkeys(element_list, 0.)
    mol_frac_element = dict.fromkeys(element_list, 0.)
    moles_g = dict.fromkeys(element_list, 0.)
   

    for i in range(len(element_list)):
        e = element_list[i]
        if (e == 'Fe'):
            continue
        elif (e == 'Si' or e == 'Ni' or e == 'S' or e == 'O' or e == 'C'):
            solutetotal = solutetotal + metalconc[e]
            mol_prop_element[e] = metalconc[e]/atdata[e]
            total = total + mol_prop_element[e]
        else:
            solutetotal = solutetotal + metalconc[e]*1e-4
            mol_prop_element[e] = (metalconc[e]*1e-4)/atdata[e]
            total = total + mol_prop_element[e]
        
   
                
    metalconc['Fe'] = 100.-solutetotal
    
    mol_prop_element['Fe'] = metalconc['Fe']/atdata['Fe']
    total = total + mol_prop_element['Fe']
  
    for i in range(len(element_list)):
        e = element_list[i]
        mol_frac_element[e]= mol_prop_element[e]/total

    mol_frac_element['Fe']= mol_prop_element['Fe']/total
    
    # moles_g no -  error attached, just nominal value.  Good for checking
    for i in range(len(element_list)):
        e = element_list[i]
        try:
            moles_g[e] = (unp.nominal_values(mol_frac_element[e])).item(0)
        except:
            moles_g[e] = mol_frac_element[e]
    
   

    return mol_frac_element, metalconc, moles_g

def els_to_oxides(wts, at_el_to_oxide):
    ox_wt = dict.fromkeys(wts, 0.)
    
    els = list(wts.keys())
    
    for i in els:
        ox_wt[i] = wts[i]/at_el_to_oxide[i]
        
    return ox_wt

def oxides_to_els(oxides, at_el_to_oxide):
    weights = dict.fromkeys(oxides, 0.)
    for key in oxides:
        
        weights[key]  = oxides[key] * at_el_to_oxide[key]
        
    return weights
    
    

def normalise_(oxides):
    
    f1 = sum(oxides.values()) - oxides['Fe']
    f2 = 100.-oxides['Fe']
    
    for key in oxides:
        if key == 'Fe':
            continue
        oxides[key] *=  (f2/f1)
    
    
    return oxides
        