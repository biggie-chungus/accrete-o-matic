#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:47:08 2019

@author: jonw
"""

def first_step(step, cf, ms_D, element_list, metalconc, f_ol, f_liq, ol_D, liq):
    """
    Calculate element distribution after first step.
    
    At each step, equilibrate mantle liquid with infalling material of size [step].
    Calculate new olivine composition and add to previous, and new liquid composition and add to previous (mant[]).
    """
    # Ensure metalconc is a dictionary
    if not isinstance(metalconc, dict):
        raise TypeError("metalconc must be a dictionary")

    sil_conc = {element: 0.0 for element in element_list}
    core = {element: 0.0 for element in element_list}
    mant = {element: 0.0 for element in element_list}
    ave_d = {element: 0.0 for element in element_list}
    ol = {element: 0.0 for element in element_list}

    for element in element_list:
        if element in {'Fe', 'O', 'S', 'C'}:
            continue
        
        # Calculate liquid concentration
        liq[element] = (liq[element] * step * f_liq + step) / (step * f_ol * ol_D[element] + step * f_liq + step * cf * ms_D[element])
        
        # Calculate olivine concentration
        ol[element] = (step - f_liq * liq[element] * step - cf * ms_D[element] * liq[element] * step) / (f_ol * step)
        
        # Calculate metal concentration
        metalconc[element] = ms_D[element] * liq[element]
        
        # Instantaneous mantle concentration
        sil_conc[element] = (f_ol * ol[element] + f_liq * liq[element]) / (f_ol + f_liq)
        
        # Update core and mantle concentrations
        core[element] += step * metalconc[element]
        mant[element] += step * sil_conc[element]
        
        # Calculate average distribution coefficient
        ave_d[element] = core[element] / mant[element]
        
    return metalconc, sil_conc, core, mant, ave_d, liq, ol