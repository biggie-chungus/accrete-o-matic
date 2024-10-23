import numpy as np
import pandas as pd
import csv
import uncertainties.unumpy as unp
import uncertainties as unc
from uncertainties import ufloat
import numpy.polynomial.polynomial as poly
from scipy.optimize import root

from typing import Dict, List, Tuple
from dataclasses import dataclass

import math
import copy
import time

from first_step import first_step

from molefrac import molefrac_
from molefrac import els_to_oxides
from molefrac import oxides_to_els
from molefrac import normalise_

from Olivine_ import Olivine_ 
from Olivine_ import Silicate_moles


"""
@ lines 717 to 713 chnage pressures, amount of Olivine C & S contents of metal
"""

@dataclass
class GeochemParams:
    """Store geochemical parameters."""
    pressure: float
    temperature: float
    oxygen_fugacity: float

def e_to_eps(gdata: Dict[str, float], atdata: Dict[str, float]) -> Dict[str, float]:
    """Convert interactions to epsilon values."""
    eps_values = {}
    for key, value in gdata.items():
        atmass = atdata[key]
        eps = ((((value/0.004343) - 1.) * atmass) / 55.85) + 1
        eps_values[key] = eps
    return eps_values

def eps_temp(temp: float, ep1875_df: pd.DataFrame, element_list: List[str]) -> pd.DataFrame:
    """Correct epsilons for temperature."""
    eps_df = ep1875_df.copy()
    eps_df = eps_df.apply(np.sqrt)
    return eps_df

def eps_elements_temp(eps_df, temp):
    """
    Returns the epsilon list, temperature corrected.
    """
    # Deep copy of the dataframe to avoid modifying the original
    epsilons = copy.deepcopy(eps_df)
    
    # Iterate over the dataframe indices
    for ei in epsilons.index:
        for ej in epsilons.index:
            # Update the epsilon values with temperature correction
            epsilons.loc[ei, ej + '_j'] = (1873.0 / temp) * eps_df.loc[ei, ej + '_j']

    return epsilons


def gam_temp(temp, atdata, element_list, gdata):
    """
    correct each gamma0 to experimental temperature
    """

    glog = 0
    gm0 = []
    for els in range(len(element_list)):

        if (element_list[els] == 'Fe'):
                continue
        else:
            mass = element_list[els]
            glog = 1873./temp*np.log(gdata[mass])
            g0 = np.exp(glog)

            gm0.append(g0)
    return gm0




def gammasolvent(temp, moles, element_list, temp_eps_df):
    """
    Calculate the gamma solvent value.
    """
    gams = 0.0
    el_number = len(element_list)

    for xi in range(el_number):
        eli = element_list[xi]
        moles[eli] = moles[eli].nominal_value if hasattr(moles[eli], 'nominal_value') else moles[eli]

        for xj in range(el_number):
            if xi == xj:
                tc_eii = temp_eps_df.loc[eli, eli + '_j']
                gams += tc_eii * (moles[eli] + math.log(1.0 - moles[eli]))
                continue

    for jm in range(el_number - 1):
        for kn in range(jm + 1, el_number):
            xj = element_list[jm]
            xk = element_list[kn]
            tc_jk = temp_eps_df.loc[xj, xk + '_j']

            if tc_jk == 0:
                continue

            gams -= tc_jk * moles[xj] * moles[xk] * (1.0 + math.log(1.0 - moles[xj]) / moles[xj] + math.log(1.0 - moles[xk]) / moles[xk])
            gams += 0.5 * tc_jk * moles[xj] * moles[xj] * moles[xk] * moles[xk] * (1.0 / (1.0 - moles[xj]) + 1.0 / (1.0 - moles[xk]) - 1.0)

    for ji in range(el_number):
        for jk in range(el_number):
            yi = element_list[ji]
            yk = element_list[jk]

            if yi == yk:
                continue

            tc_jk = temp_eps_df.loc[yi, yk + '_j']

            if tc_jk == 0:
                continue

            gams += tc_jk * moles[yi] * moles[yk] * (1.0 + math.log(1.0 - moles[yk]) / moles[yk] - 1.0 / (1.0 - moles[yi]))
            gams -= tc_jk * moles[yi] * moles[yi] * moles[yk] * moles[yk] * (1.0 / (1.0 - moles[yi]) + 1.0 / (1.0 - moles[yk]) + moles[yi] / (2.0 * (1.0 - moles[yi]) ** 2) - 1.0)

    return math.exp(gams)

def gammasolute(moles, gamma_solvent, gammadata, element_list, temp_eps_df):
    """
    Calculate the gamma solute values.
    """
    gamma_solutes = {}

    for ei in element_list:
        solute_e = 0.0
        solute_f = 0.0
        solute_g = 0.0

        tcsolute_eii = temp_eps_df.loc[ei, ei + '_j']
        solute_e = math.log(gamma_solvent) + math.log(gammadata[ei]) - tcsolute_eii * math.log(1.0 - moles[ei])

        for fk in element_list:
            if ei == fk:
                continue

            if fk in {'Si', 'Ni', 'S', 'O', 'C'}:
                solute_f += temp_eps_df.loc[ei, fk + '_j'] * moles[fk] * (1.0 + math.log(1.0 - moles[fk]) / moles[fk] - 1.0 / (1.0 - moles[ei]))

        for gk in element_list:
            if ei == gk:
                continue

            if gk in {'Si', 'Ni', 'S', 'O', 'C'}:
                solute_g += temp_eps_df.loc[ei, gk + '_j'] * moles[gk] ** 2 * moles[ei] * (1.0 / (1.0 - moles[ei]) + 1.0 / (1.0 - moles[gk]) + moles[ei] / (2.0 * (1.0 - moles[ei]) ** 2) - 1.0)

        solute_gamma = solute_e - solute_f + solute_g
        gamma_solutes[ei] = math.exp(solute_gamma)

    return gamma_solutes






def carbon_content(x, temp):
    """
    calculate metal carbon saturation from Wood 1993 and references therein
    optimise this for C content at C-saturation.
    Go on, I dare you.  Do it.

    """
    if temp > 2250.:
        g0cc = -14885. - 64.66*(temp-2250.)
    else:
        g0cc = 22600.-16.66*temp

    thetacc = 7830./temp+1.658


    fit = g0cc + (8.314*temp*np.log(x/(1.-2.*x))
                 + 8.314*temp*(thetacc*x/(1.-x)))
    return fit




def carbon_calculator(atdata, moles, temp):
    """
    Calculates the carbon content of the metal at C-saturation.
    Uses data from Wang et al 1991 and Wood 1993.
    """
    def carbon_content_func(x):
        return carbon_content(x, temp)

    # Use root to find the solution
    solution = root(carbon_content_func, 0.2)
    res = solution.x[0]

    percent_c = (res * atdata['C'] / (((1.0 - res) * atdata['Fe']) + res * atdata['C'])) * 100

    # Adjust Fe moles
    adjusted_Fe = (1.0 - res) * initial_moles_fe
    mult = adjusted_Fe / initial_moles_fe

    # Update moles for each element
    for element in moles:
        moles[element] *= mult
    moles['C'] = res * initial_moles_fe
    moles['Fe'] = adjusted_Fe

    return moles, percent_c




def calculate_nbot(oxide_moles):
    Ca_O = oxide_moles['Ca']
    Mg_O = oxide_moles['Mg']
    Fe_O = oxide_moles['Fe']
    Na2_O = oxide_moles['Na']
    K2_O = oxide_moles['K']
    Si_O2 = oxide_moles['Si']
    Ti_O2 = oxide_moles['Ti']
    Al2_O3 = oxide_moles['Al']

    yNBOT = 2 * (Ca_O + Mg_O + Fe_O + Na2_O + K2_O - Al2_O3)
    xT = Si_O2 + 2 * Al2_O3 + Ti_O2
    return yNBOT / xT

def calculate_altnbot(element_moles):
    t = element_moles['Al'] + element_moles['Si']
    o = 1.0 - sum(element_moles.values())
    altnbot = (2 * o - 4 * t) / t
    return max(0, min(altnbot, 4))

def calculate_ms_D(e, gamma_solvent, gamma_solutes, moles, oxide_moles, element_moles, metal_fe, silicate_ox, 
                   temp, pressure, si_melt, fe_melt, ni_melt, w_melt, altnbot):
    """
    

    Parameters
    ----------
    e : TYPE
        DESCRIPTION.
    gamma_solvent : TYPE
        DESCRIPTION.
    gamma_solutes : TYPE
        DESCRIPTION.
    moles : TYPE
        DESCRIPTION.
    oxide_moles : TYPE
        DESCRIPTION.
    element_moles : TYPE
        DESCRIPTION.
    temp : TYPE
        DESCRIPTION.
    pressure : TYPE
        DESCRIPTION.
    si_melt : TYPE
        DESCRIPTION.
    fe_melt : TYPE
        DESCRIPTION.
    ni_melt : TYPE
        DESCRIPTION.
    w_melt : TYPE
        DESCRIPTION.
    altnbot : TYPE
        DESCRIPTION.

    Returns
    -------
    ms_D : TYPE
        DESCRIPTION.

    """
    
    ms_D = 0
    gamma_solvent_fe = gamma_solvent * unp.nominal_values(moles['Fe']) / unp.nominal_values(oxide_moles['Fe'])
    
   
    
    # gamma_solvent_fe_w = gamma_solvent * unp.nominal_values(moles['Fe']) / unp.nominal_values(element_moles['Fe'])
    if e == 'Si':
        b = si_melt
        gsi = math.log10(gamma_solutes['Si'])
        ms_D = 10 ** (b / temp - gsi + 2 * math.log10(gamma_solvent_fe))
    
    elif e == 'Co':
        gco = math.log10(gamma_solutes['Co'])
        if pressure <= 6:
            ms_D = 10 ** (ufloat(1.287, 0.163) + ufloat(-227., 288.) / temp - gco + ufloat(-170.265, 25.712) * pressure / temp + math.log10(gamma_solvent_fe))
        else:
            ms_D = 10 ** (ufloat(0.3, 0.18) + ufloat(1405., 369.) / temp - gco + ufloat(-40., 7.) * pressure / temp + math.log10(gamma_solvent_fe))
    
    elif e == 'W':
        gw = math.log10(gamma_solutes['W'])
        # ms_D = 10 ** (ufloat(1.85, 0.24) + (w_melt) / temp + ufloat(-77., 10) * pressure / temp - gw + 3 * math.log10(gamma_solvent_fe))
        ms_D = 10 ** (ufloat(1.85, 0.24) + (w_melt+2637) / temp + ufloat(-77., 10) * pressure / temp - gw + 3 * math.log10(gamma_solvent_fe))
        
        
    elif e == 'Ni':
        gni = math.log10(gamma_solutes['Ni'])
        if pressure <= 6:
            ms_D = 10 ** (ufloat(0.001, 0.220) + ufloat(3503, 382.) / temp - gni - math.log10(fe_melt / ni_melt) + ufloat(-185.606, 36.896) * pressure / temp + math.log10(gamma_solvent_fe))
    
    elif e == 'Mo':
        gmo = math.log10(gamma_solutes['Mo'])
        ms_D = 10 ** (ufloat(0.42, 0.1) + ufloat(3117, 539) / temp - gmo + ufloat(-97, 25.) * pressure / temp + ufloat(-0.27, 0.03) * altnbot + 2 * math.log10(gamma_solvent_fe))
    
    elif e == 'Cr':
        """ 
        nb wood et al. use wt D's not MF D's ,Tuff et al 2010 wood 2008 eqn13
        """
        gcr = math.log10(gamma_solutes['Cr'])
        ms_D = 10 ** (ufloat(0.643, 0.1) + ufloat(-4323, 538) / temp - gcr + ufloat(-22, 13) * pressure / temp + math.log10(gamma_solvent_fe))
    
    elif e == 'Mn':
        gmn = math.log10(gamma_solutes['Mn'])
        ms_D = 10 ** (ufloat(0.04, 0.1) + ufloat(-5761, 576) / temp - gmn + ufloat(-49, 16) * pressure / temp + math.log10(gamma_solvent_fe))
    
    elif e == 'Nb':
        """ 
        nb wood et al. use wt D's not MF D's
        """
        gnb = math.log10(gamma_solutes['Nb'])
        ms_D = 10 ** (ufloat(2.837, 0.41) + ufloat(-15500, 2000) / temp - gnb + ufloat(-114., 43) * pressure / temp + ufloat(-0.47, .30) * altnbot + 2.5 * math.log10(gamma_solvent * metal_fe/ silicate_ox['Fe']))
    
    elif e == 'V':
        """ 
        nb wood et al. use wt D's not MF D's
        """
        gv = math.log10(gamma_solutes['V'])
        ms_D = 10 ** (ufloat(0.855, 0.1) + ufloat(-8548, 854) / temp - gv + ufloat(-62, 19) * pressure / temp + ufloat(-0.101, 0.029) * altnbot + 1.5 * math.log10(gamma_solvent * (metal_fe / silicate_ox['Fe'])))
        
    return ms_D

def element_d(element_list, gamma_solvent, gamma_solutes, oxide_moles, element_moles, temp, pressure, moles,  metalconc, silicate_ox, w_melt, si_melt, fe_melt, ni_melt):
    ms_D = {element: 0.0 for element in element_list}
    
    try:
        metal_fe = metalconc['Fe'].nominal_value
    except AttributeError:
        metal_fe = metalconc['Fe']
    
    nbot = calculate_nbot(oxide_moles)
    altnbot = calculate_altnbot(element_moles)
    
    for element in element_list:
        ms_D[element] = calculate_ms_D(element, gamma_solvent, gamma_solutes, moles, oxide_moles, element_moles, metal_fe, silicate_ox,temp, pressure, si_melt, fe_melt, ni_melt, w_melt, altnbot)
    
    return ms_D, altnbot

def calculate_volume(temp, pressure, v298, k0, kprime, a, b, c, d):
    """
    Calculate the volume at given temperature and pressure.
    using the Ozawa et al model for O partitioning between fp & metal
    """
    vt = (1000.0 * v298 * np.exp(a * (temp - 298.0) + 0.5 * b * (temp**2 - 298.0**2) + c * math.log(temp / 298.0) + d / temp - d / 298.0))
    volume = vt * (k0 / (kprime - 1.0) * ((1.0 + kprime * pressure / k0)**((kprime - 1.0) / kprime) - 1.0))
    return volume

def o_(i, step, pressure, temp, ave_d, oxide_moles, core, mant, metalconc, moles, sflag, gOratio, atdata):
    """
    not a big deal at low pressures.
    mfeo = the molefraction of FeO in silicate, corrected for partitioning
    between mw and silicate melt - see Dave Rubie's paper
    This is a bit mucky but is here because its used to calculate the oxygen
    content in metal during an accretion run
    """
    
    mfeo = (1.148 * oxide_moles['Fe'] * 0.0094 + 1.319 * (oxide_moles['Fe'] * 0.0094)**2)

    # using the Ozawa et al model for O partitioning between fp & metal
    """
    magnesiowustite  - labelled solid in paper
    """
    Vmw = calculate_volume(temp, pressure, 12.25, 180.0, 4.0, 2.6e-5, 1.47e-8, 2.79e-3, 4.28)
    
    """
    oxygen in liquid Fe
    """
    Vox = calculate_volume(temp, pressure, 5.99, 103.7, 4.0, 1.06e-5, 7.16e-10, 2.79e-5, 0.0)
    
    """
    Liquid Fe
    """
    Vfe = calculate_volume(temp, pressure, 6.92, 109.7, 4.661, 9.27e-5, 0.0, 0.0, 0.0)
    
    """
    above checked  - and seems to work according to paper
    """
    dH = 170000.0
    dS = 55.0

    lnkd = (-dH + temp * dS - (Vfe + Vox - Vmw)) / (8.31444621 * temp)
    kd = np.exp(lnkd)
    
    """
    sflag == 1, sulfur is present in the metal and g0ratio applied
    if not, don't apply it.  Obvs.
    """
    if sflag > 0:
        o_coremf = kd * (mfeo / moles['Fe']) * gOratio
    else:
        o_coremf = kd * (mfeo / moles['Fe'])
    
    """
    convert mole fraction to wt% for oxygen
    moles of Fe kept constant since paper only considers partitioning btwn
    Fe and FeO
    """
    numo = o_coremf * 15.9994
    numfe = moles['Fe'] * 55.845
    O_percent = (numo / (numo + numfe) * 100)

    """
    the upper limit of oxygen in core is set by the FeO
    content of the silicate- if calculated FeO in metal
    is greater than FeO content of silicate, then metallic FeO content
    is forced to be the same as the silicate FeO content
    """
    if O_percent > (oxide_moles['Fe'] * 1.2865):
        metalconc['O'] = oxide_moles['Fe']
    else:
        metalconc['O'] = O_percent

    core['O'] += step * metalconc['O']
    mant['O'] += step * mfeo

    if i > 0:
        ave_d['O'] = core['O'] / mant['O']
    else:
        ave_d['O'] = 0

    return metalconc, core, mant, ave_d

def calculate_oxygen_fugacity(temp, pressure):
    """
    Calculate the oxygen fugacity.
    """
    Pbar = pressure * 10000
    nno = (-24930 / temp) + 9.36 + 0.046 * (Pbar - 1) / temp
    iw = (-27489 / temp) + 6.702 + 0.055 * (Pbar - 1) / temp
    return nno, iw

def calculate_effective_oxygen_fugacity(oxide_moles, gamma_solvent, moles):
    """
    Calculate the  oxygen fugacity relative to IW.
    """
    effO2 = 2 * math.log10(oxide_moles['Fe'] / (gamma_solvent * unp.nominal_values(moles['Fe'])))
    return effO2

def calculate_partition_coefficients(diw, dnno):
    """
    Calculate the partition coefficients for V and Cr betweeen silicacte melt and olivine
    taken from Malmann & O'neill 2009
    """
    ol_D = {}
    ol_D['V'] = 10 ** ((-0.24) * dnno - 1.41)
    ol_D['Cr'] = (8.5e-1 * 5.95e-3 ** -1 * (10 ** diw) ** -0.5) / ((5.95e-3 ** -1 * (10 ** diw) ** -0.5) + 1)
    return ol_D

def next_step(step, count, f_met, ms_D, element_list, core, ave_d, f_ol, ol_D, f_liq, pressure, temp, oxide_moles, gamma_solvent, moles, liq, metalconc, mant):
    """
    Calculate element distribution after first step.
    
    Calculates D(ol/melt) for elements that are dependent upon fO2.
    
    At each step, equilibrate mantle liquid with infalling material of size [step].
    Calculate new olivine composition and add to previous, and new liquid composition and add to previous (mant[]).
    """
    ol = {element: 0.0 for element in element_list}
    sil_conc = {element: 0.0 for element in element_list}
    
    nno, iw = calculate_oxygen_fugacity(temp, pressure)
    effO2 = calculate_effective_oxygen_fugacity(oxide_moles, gamma_solvent, moles)
    diw = iw + effO2
    dnno = diw - nno
    
    ol_D.update(calculate_partition_coefficients(diw, dnno))
    
    for element in element_list:
        if element in {'Fe', 'O', 'S', 'C'}:
            continue
        
        liq[element] = (liq[element] * count * f_liq + step) / (step * f_ol * ol_D[element] + (count + step) * f_liq + step * f_met * ms_D[element])
        ol[element] = (step - f_liq * liq[element] * step - f_met * ms_D[element] * liq[element] * step) / (f_ol * step)
        metalconc[element] = ms_D[element] * liq[element]
        
        # Instantaneous mantle concentration
        sil_conc[element] = (f_ol * ol[element] + f_liq * liq[element]) / (f_ol + f_liq)
        core[element] += step * metalconc[element]
        mant[element] += step * sil_conc[element]
        
        ave_d[element] = core[element] / mant[element]
        
    return metalconc, sil_conc, core, mant, ave_d, liq, ol

def Silicate_elements(silicate_els, silicate, ave_d, corefraction, metalconc, chondrite):
    """
    Returns the elemental mantle concentration relative to chondrite for major elements.
    """
    for element in silicate_els:
        if element in {'Fe', 'Si'}:
            continue
        
        silicate[element] = 1.0 / (corefraction + (1 - corefraction)) * chondrite[element]

    return silicate


def calculate_activity(oxide_moles, coefficients, t):
    """
    Calculate the activity of a component in the melt.
    """
    
    # Sum all coefficients
    activity = sum(coefficients.values())
    return np.exp(activity / (8.314 * t))

def ni_melt_(oxide_moles, t):
    """
    Calculate the activity of Ni in the melt from wood 2013
    """
    coefficients = {
        'const': 12332.,
        'Si_Mg': -55097. * oxide_moles['Si'] * oxide_moles['Mg'],
        'Ca_Mg': 91837 * oxide_moles['Ca'] * oxide_moles['Mg'],
        'Na_K_Si': 299721 * (oxide_moles['Na'] + oxide_moles['K']) * oxide_moles['Si'],
        'Na_K_Al': -580839 * (oxide_moles['Na'] + oxide_moles['K']) * oxide_moles['Al'],
        'Ti': -60241 * oxide_moles['Ti']
    }
    
    return calculate_activity(oxide_moles, coefficients, t)

def fe_melt_(oxide_moles, t):
    """
    Calculate the activity of Fe in the melt
    """
    coefficients = {
        'const': -18994.,
        'Na': 99579. * oxide_moles['Na'],
        'Al_Si': 192553 * oxide_moles['Al'] * oxide_moles['Si'],
        'Ca_Mg': 282789. * oxide_moles['Ca'] * oxide_moles['Mg'],
        'Ca': 79492 * oxide_moles['Ca'] ** 2,
        'Fe': 120972 * oxide_moles['Fe'] ** 2
    }
    return calculate_activity(oxide_moles, coefficients, t)

def w_melt_(oxide_moles, t):
    const  = unc.ufloat(-73730, 6030)
    si = unc.ufloat(7.4e4, 1.2e4) * oxide_moles['Si']**2
    mg = unc.ufloat(-6.3e4, 1.2e4) * oxide_moles['Mg']**2
    mgsi = unc.ufloat(1.1e5, 0.3e5) * oxide_moles['Si']*oxide_moles['Mg']
    ca = unc.ufloat(-7.3e5, 1.0e5) * oxide_moles['Ca']**2
    fe = unc.ufloat(1.73e5, 0.49e5) * oxide_moles['Fe']**2
    alca = unc.ufloat(2.2e6, 0.3e6) * oxide_moles['Al']*oxide_moles['Ca']
    melt = (const + si + mg + mgsi + ca + fe + alca)/8.314
    
    return melt


def si_melt_(oxide_moles, t):
    """
    Calculate the activity of Si in the melt.
    """
    e0 =  unc.ufloat(-77970, 4585)
    e1 = unc.ufloat(-1.824e5, 34910) * oxide_moles['Ca']*oxide_moles['Mg']
    e2 = unc.ufloat(-1.334e5, 13606) * oxide_moles['Mg']**2
    e3 = unc.ufloat(-4.202e5, 66318)  * oxide_moles['Ca']*oxide_moles['Si']
    
    return  (e0+e1+e2+e3)/8.314
    



def s_absent (moles, metalconc, atdata, temp_eps_df, element_list,
              gammadata, gamma_element_list, ms_D, i, temp, cpres):
    
    """
    used at high S contents (approaching FeS) where  activity of Fe is poorly constrained
    this and s_present just places 'reasonable' constraints on the activity of Fe
    This isn't an issue here where S content is fixed and low (<say, 6wt%)
    """

    metalconc_sabs = dict.fromkeys(element_list, 0.000001)
    metalconc_sabs['Si'] = metalconc['Si']
    metalconc_sabs['O'] = metalconc['O']
    metalconc_sabs['C'] = metalconc['C']
    metalconc_sabs['Fe'] = metalconc['Fe'].nominal_value
    

    sabs_moles, sabs_metalconc, moles_g_abs = molefrac_(atdata, element_list,  metalconc_sabs, cpres)
    gamma_solvent = gammasolvent(temp, sabs_moles, gamma_element_list,
                                 temp_eps_df)
    
    gammacheck = gamma_solvent * moles['Fe']
    if (i < 1. and gammacheck > 1.):
        
        gamma_solvent = 1./moles['Fe']
        print ('SeaSponge!!',moles['Fe'])
       
    gamma_solutes = gammasolute(sabs_moles, gamma_solvent, gammadata,
                                       gamma_element_list, temp_eps_df)
    gOa = gamma_solutes['O']

    return gOa

def s_present (moles, metalconc, atdata, temp_eps_df, element_list,
              gammadata, gamma_element_list, ms_D, i, temp, cpres):
    
    """
    used at high S contents (approaching FeS) where  activity of Fe is poorly constrained
    see above
    """

    metalconc_spres = dict.fromkeys(element_list, 0.000001)
    metalconc_spres['Si'] = metalconc['Si']
    metalconc_spres['O'] = metalconc['O']
    metalconc_spres['C'] = metalconc['C']
    metalconc_spres['S'] = metalconc['S']

    spres_moles, spres_metalconc, moles_g_pres = molefrac_(atdata, element_list,  metalconc_spres, cpres)
    gamma_solvent = gammasolvent(temp, spres_moles, gamma_element_list,
                                 temp_eps_df)

    gammacheck = gamma_solvent * moles['Fe']
    if (i < 1. and gammacheck > 1.):
        gamma_solvent = 1./moles['Fe']
#        print('SECOND_C0CK! (line 646)')
    gamma_solutes = gammasolute(spres_moles, gamma_solvent, gammadata,
                                       gamma_element_list, temp_eps_df)
    gOs = gamma_solutes['O']

    return gOs




def read_csv_to_dict(file_path, key_col, value_col):
    data = np.genfromtxt(file_path, delimiter=',', dtype=str)
    keys = data[:, key_col]
    values = data[:, value_col].astype(float)
    return {key: value for key, value in zip(keys, values)}

def initialize_variables():
    fraction = []
    av_ni, av_co, av_cr, av_v, av_w, av_si, av_nb, av_mo = ([] for _ in range(8))
    inst_feo, inst_cf = [], []
    ni_melt, fe_melt, w_melt, si_melt = 0, 0, 0, 0
    return fraction, av_ni, av_co, av_cr, av_v, av_w, av_si, av_nb, av_mo, inst_feo, inst_cf, ni_melt, fe_melt, w_melt, si_melt

def read_element_data():
    atom_elems = np.genfromtxt('atmass.csv', delimiter=',', dtype=str)[:, 0]
    atmass = np.genfromtxt('atmass.csv', delimiter=',')[:, 1]
    atdata = {elem: mass for elem, mass in zip(atom_elems, atmass)}

    atmulti = np.genfromtxt('atmass.csv', delimiter=',')[:, 2]
    at_el_to_oxide = {elem: multi for elem, multi in zip(atom_elems, atmulti)}

    v = np.genfromtxt('atmass.csv', delimiter=',')[:, 4]
    valence = {elem: val for elem, val in zip(atom_elems, v)}

    return atdata, at_el_to_oxide, valence

def read_gamma_data():
    elements = np.genfromtxt('gamma0.csv', delimiter=',', usecols=0, dtype=str)
    gamma0s = np.genfromtxt('gamma0.csv', delimiter=',')[:, 1]
    gdata = {elem: gamma for elem, gamma in zip(elements, gamma0s)}
    return gdata

def read_epsilon_data():
    return pd.read_csv('epsilons.csv', index_col=0, na_values=['nan', ''])

def initialize_output_file():
    out = open('accretion_oxidised_1.csv', 'w')
    head = csv.writer(out, delimiter=',')
    head.writerow(['Pressure', 'Temp', 'nbot', 'fe_mant', 'feO_mant', 'S_metal', 'silicate_olivine%', 
                   'ol_frac', 'C_metal', 'cf', 'Liq_frac', 'ave_ni', 'ni_err', 'ave_co', 'co_err', 
                   'ave_cr', 'cr_err', 'ave_nb', 'nb_err', 'ave_w', 'w_err', 'ave_v', 'v_err', 
                   'ave_mo', 'mo_err', 'ave_si', 'si_err', 'ave_mn', 'mn_err', 'ni_act', 'fe_act', 
                   'Ol_D_Cr', 'Ol_D_V', 'Ol_D_Ni', 'Ol_D_Co'])
    return out

def main():
    fraction, av_ni, av_co, av_cr, av_v, av_w, av_si, av_nb, av_mo, inst_feo, inst_cf, ni_melt, fe_melt, w_melt, si_melt = initialize_variables()
    atdata, at_el_to_oxide, valence = read_element_data()
    gdata = read_gamma_data()
    eps_df = read_epsilon_data()

    
    gamma_element_list = list(element_list)
    gamma_element_list.remove('Fe')

    chondrite = {'Si': 19.2, 'Ti': 0.137, 'Al': 1.566, 'Mg': 17.37, 'Fe': 32.6, 'Ca': 1.7, 'Na': .178, 'K': 0.016}
   

    temp_melts = [1591.95, 1589.61, 1581.99, 1571.45, 1552.70]
    feo_melts = [2., 5., 10., 15., 22.]
    temp_k = [x + 273 for x in temp_melts]
    temp_coeff = poly.polyfit(feo_melts, temp_k, 3)

    

    out = initialize_output_file()
    # n.b. split the pressures and merge the output to speed things up......
    f_ = [ 10,12,14,16, 18, 20, 22, 24, 26, 28, 30]                  # mantle Fe content (wt%)
    p_splat = [3]                                                    # Final pressure of segregation - can loop over list e.g. [0.5, 1, 3 etc...]
    ol_ = [1, 10, 20, 30, 40, 50, 60]                                # percent of melt that is Olivine 

    c_= [1e-3, 1, 2, 3, 4, 5, 6]                                     # C in metal (wt%), 6 is assumed C-saturated
    s_ = [1e-3, 1, 2, 3, 4, 5, 6]    					     		  # S in metal (wt%)

    for peak_pressure in p_splat:
        for feend in f_:
            
            for slph in s_:
                for ol_percent in ol_:
                    for c in c_:
                        i, iold, step, counter, increment_flag = 0.0, 0.0, 0.01, 0, 0
                        pp, total, total_cf, total_feo = 1.0 - step, 0.0, 0.0, 0.0

                        silicate_els = ['Si', 'Ti', 'Al', 'Mg', 'Fe', 'Ca', 'Na', 'K']
                        sil_wts = dict.fromkeys(silicate_els, 0.0)
                        sil = [21., 0.12, 2.35, 22.8, 6.26, 2.53, 0.26, 0.024]
                        sil_wts_bulk = {elem: wt for elem, wt in zip(silicate_els, sil)}
                        sil_ox_bulk = els_to_oxides(sil_wts_bulk, at_el_to_oxide)
                        oxide_moles = dict.fromkeys(silicate_els, 0.0)

                        ave_d = dict.fromkeys(element_list, 0.1)
                        metalconc = dict.fromkeys(element_list, 0.001)
                        sil_conc = dict.fromkeys(element_list, 0.0)
                        liq = dict.fromkeys(element_list, 0.0)
                        elements = dict.fromkeys(element_list, 0.01)
                        core = dict.fromkeys(element_list, 0.0)
                        mant = dict.fromkeys(element_list, 0.0)
                        moles = dict.fromkeys(element_list, 0.0)

                        

                        sil_wts_bulk['Fe'] = feend
                        sil_ox_bulk = els_to_oxides(sil_wts_bulk, at_el_to_oxide)
                        sil_ox_bulk = normalise_(sil_ox_bulk)
                        sil_wts_bulk = oxides_to_els(sil_ox_bulk, at_el_to_oxide)

                        temp = (temp_coeff[0] + sil_ox_bulk['Fe'] * temp_coeff[1] +
                                sil_ox_bulk['Fe'] ** 2 * temp_coeff[2] + sil_ox_bulk['Fe'] ** 3 * temp_coeff[3]) + 28.57 * 0.01 - ol_percent

                        if c == 6:
                            cpres = 'y'
                            moles, percent_c = carbon_calculator(atdata, moles, temp)
                            metalconc['C'] = percent_c
                            
                        else:
                            cpres = 'n'
                            metalconc['C'] = c
                            

                        spres = 'y'
                        sflag = 1
                        metalconc['S'] = slph

                        moles, metalconc, moles_g = molefrac_(atdata, element_list, metalconc, cpres)

                        cf = (31.5 - sil_wts_bulk['Fe']) / (85 - sil_wts_bulk['Fe'])

                        while i < (1 + step):
                            if i < step:
                                print('Max pressure = ', peak_pressure, 'GPa')
                                print('FeO start', sil_ox_bulk['Fe'], '%')
                                print('FeO end', feend / at_el_to_oxide['Fe'], '%')
                                print('olivine % = ', ol_percent)
                                print('carbon % = ', c)
                                print('sulphur % = ', metalconc['S'])

                                total_cf += step * cf
                                total_feo += step * sil_wts_bulk['Fe']
                                total += step

                                silicate_frac = 1.0 - cf
                                f_ol = (silicate_frac * ol_percent) / 100.0
                                f_liq = silicate_frac - f_ol

                                pressure = peak_pressure * ((i / pp) ** 0.2) ** 2

                                temp_eps_df = eps_elements_temp(eps_df, temp)
                                gamma0 = gam_temp(temp, atdata, gamma_element_list, gdata)
                                gammadata = {elem: gamma for elem, gamma in zip(gamma_element_list, gamma0)}

                                gamma_solvent = gammasolvent(temp, moles, gamma_element_list, temp_eps_df)
                                gamma_solutes = gammasolute(moles_g, gamma_solvent, gammadata, gamma_element_list, temp_eps_df)

                                sil_ox, ol_D, olmg = Olivine_(ol_percent, sil_ox_bulk, element_list)
                                sil_wts = oxides_to_els(sil_ox, at_el_to_oxide)
                                ox_sum = sum(sil_ox_bulk.values())

                                element_moles, oxide_moles = Silicate_moles(sil_wts, ox_sum, silicate_els, at_el_to_oxide, atdata, valence)

                                w_melt = w_melt_(oxide_moles, temp)
                                si_melt = si_melt_(oxide_moles, temp)
                                ni_melt = ni_melt_(oxide_moles, temp)
                                fe_melt = fe_melt_(oxide_moles, temp)

                                ms_D, altnbot = element_d(element_list, gamma_solvent, gamma_solutes, oxide_moles, element_moles, temp, pressure,  moles, metalconc, sil_ox, w_melt, si_melt, fe_melt, ni_melt)

                                metalconc, core, mant, ave_d = o_(i, step, pressure, temp, ave_d, oxide_moles, core, mant, metalconc,  moles, sflag, 1, atdata)
                                
                                metalconc, sil_conc, core, mant, ave_d, liq, ol = first_step(step, cf, ms_D, element_list, metalconc, f_ol, f_liq, ol_D, liq)

                                moles, metalconc, moles_g = molefrac_(atdata, element_list, metalconc, cpres)
                                silicate = Silicate_elements(silicate_els, sil_wts, ave_d, cf, metalconc, chondrite)

                            if i > step:
                                fe_silicate = feend
                                sil_wts_bulk['Fe'] = fe_silicate
                                sil_ox_bulk = els_to_oxides(sil_wts_bulk, at_el_to_oxide)
                                sil_ox_bulk = normalise_(sil_ox_bulk)
                                sil_wts_bulk = oxides_to_els(sil_ox_bulk, at_el_to_oxide)

                                for key in sil_wts:
                                    sil_wts[key] = sil_wts[key] * (1.0 - step) + step * sil_wts_bulk[key]

                                for key in sil_ox:
                                    sil_wts[key] = sil_ox[key] * at_el_to_oxide[key]

                                cf = (31.5 - sil_wts_bulk['Fe']) / (85 - sil_wts_bulk['Fe'])
                                silicate_frac = 1.0 - cf
                                f_ol = (silicate_frac * ol_percent) / 100.0
                                f_liq = silicate_frac - f_ol

                                total_feo += step * sil_wts_bulk['Fe']
                                total_cf += step * cf
                                total += step

                                if i > 0.99:
                                    pressure = peak_pressure
                                else:
                                    pressure = peak_pressure * ((i / pp) ** 0.2) ** 2

                                temp = (temp_coeff[0] + sil_ox['Fe'] * temp_coeff[1] +
                                        sil_ox['Fe'] ** 2 * temp_coeff[2] + sil_ox['Fe'] ** 3 * temp_coeff[3]) + 28.57 * pressure - ol_percent

                                if cpres in ['Y', 'y']:
                                    moles, percent_c = carbon_calculator(atdata, moles, temp)
                                    metalconc['C'] = percent_c
                                else:
                                    metalconc['C'] = c
                                

                                moles, metalconc, moles_g = molefrac_(atdata, element_list, metalconc, cpres)
                                temp_eps_df = eps_elements_temp(eps_df, temp)
                                gamma0 = gam_temp(temp, atdata, gamma_element_list, gdata)
                                gammadata = {elem: gamma for elem, gamma in zip(gamma_element_list, gamma0)}
								
                                # if spres in ['Y', 'y']:
                                #     gOa = s_absent(moles_g, metalconc, atdata, temp_eps_df, element_list, gammadata, gamma_element_list, ms_D, i, temp, cpres)
                                #     gOs = s_present(moles_g, metalconc, atdata, temp_eps_df, element_list, gammadata, gamma_element_list, ms_D, i, temp, cpres)
                                #     gOratio = gOa / gOs
                                #     gOratio = min(gOratio, 3.0)
                                # else:
                                #     gOratio = 1.0
								
                                gOratio = 1.0
                                gamma_solvent = gammasolvent(temp, moles_g, gamma_element_list, temp_eps_df)
                                gamma_solutes = gammasolute(moles_g, gamma_solvent, gammadata, gamma_element_list, temp_eps_df)

                                sil_ox, ol_D, olmg = Olivine_(ol_percent, sil_ox_bulk, element_list)
                                sil_wts = oxides_to_els(sil_ox, at_el_to_oxide)
                                ox_sum = sum(sil_ox.values())

                                element_moles, oxide_moles = Silicate_moles(sil_wts, ox_sum, silicate_els, at_el_to_oxide, atdata, valence)
                                w_melt = w_melt_(oxide_moles, temp)
                                si_melt = si_melt_(oxide_moles, temp)
                                ni_melt = ni_melt_(oxide_moles, temp)
                                fe_melt = fe_melt_(oxide_moles, temp)

                                ms_D, altnbot = element_d(element_list, gamma_solvent, gamma_solutes, oxide_moles, element_moles, temp, pressure,  moles, metalconc, sil_ox, w_melt, si_melt, fe_melt, ni_melt)
                                metalconc, core, mant, ave_d = o_(i, step, pressure, temp, ave_d, oxide_moles, core, mant, metalconc, moles, sflag, 1, atdata)

                                
                                metalconc, sil_conc, core, mant, ave_d, liq, ol = next_step(step, i, cf, ms_D, element_list, core, ave_d, f_ol, ol_D, f_liq, pressure, temp, oxide_moles, gamma_solvent, moles, liq, metalconc, mant)
                                silicate = Silicate_elements(silicate_els, silicate, ave_d, cf, metalconc, chondrite)

                            if i >= 1:
                                break
                            else:
                                i += step

                        print('total feo', total_feo)
                        print('i', i, 'temp', temp, 'feo', sil_ox_bulk['Fe'], 'ol%', ol_percent / 100.0, 'f_ol', f_ol, 'cf', cf, 'fLiq', f_liq)
                        print('ni', ave_d['Ni'], 'co', ave_d['Co'], 'mo', ave_d['Mo'])
                        print('Mn', ave_d['Mn'], 'fe_melt', fe_melt, 'ni_melt', ni_melt)
                        print('carbon content', metalconc['C'])
                        print(' ')
                        print('W', ave_d['W'])
                        print('Next!')

                        row = (pressure, '{:4.3f}'.format(temp), '{:4.3f}'.format(altnbot), '{:4.3f}'.format(sil_wts_bulk['Fe']),
                               '{:4.3f}'.format(sil_ox_bulk['Fe']), metalconc['S'], '{:5.4f}'.format(ol_percent / 100.0),
                               '{:5.4f}'.format(f_ol), '{:4.3f}'.format(metalconc['C']), '{:6.5f}'.format(cf),
                               '{:6.5f}'.format(f_liq), unp.nominal_values(ave_d['Ni']), unp.std_devs(ave_d['Ni']),
                               unp.nominal_values(ave_d['Co']), unp.std_devs(ave_d['Co']), unp.nominal_values(ave_d['Cr']),
                               unp.std_devs(ave_d['Cr']), unp.nominal_values(ave_d['Nb']), unp.std_devs(ave_d['Nb']),
                               unp.nominal_values(ave_d['W']), unp.std_devs(ave_d['W']), unp.nominal_values(ave_d['V']),
                               unp.std_devs(ave_d['V']), unp.nominal_values(ave_d['Mo']), unp.std_devs(ave_d['Mo']),
                               unp.nominal_values(ave_d['Si']), unp.std_devs(ave_d['Si']), unp.nominal_values(ave_d['Mn']),
                               unp.std_devs(ave_d['Mn']), '{:6.5f}'.format(ni_melt), '{:6.5f}'.format(fe_melt),
                               '{:6.5f}'.format(ol_D['Cr']), 
                               '{:6.5f}'.format(ol_D['V']), 
                               '{:6.5f}'.format(ol_D['Ni']), 
                               '{:6.5f}'.format(ol_D['Co']),
                               )
                        try:
                            wr = csv.writer(out,delimiter = ',')
                            wr.writerow(row)
                            time.sleep(0.5)
                        except TimeoutError as e:
                            print(f"TimeoutError: {e}")
                        except Exception as e:
                            print(f"An error occurred: {e}")
                        
    out.close()  
                    

if __name__ == "__main__":
    element_list = ['Fe', 'C', 'O', 'Si', 'Mn', 'Cr', 'V', 'Ni', 'Co', 'Nb', 'W', 'Mo', 'S']
    initial_moles_fe = 0.95
    main()

                    