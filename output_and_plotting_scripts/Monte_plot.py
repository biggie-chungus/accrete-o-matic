# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec

import random


## ==============================================================================
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
# ==============================================================================

"""
jon's big plotting script, obviously not too use friendly
firstly, read in the csv lookup table
  
change element names at lines 34, late veneer amount at l. 36 and 
change file name for plot output

will make one plot for each olivine content of oxidised melt, if so desired
"""
element = 'mo'

lv = 0.01  # late veneer
ol_frac =[0.50] # fraction of olivine in the oxidised silicate (fraction of silicate) 


# ol_frac = [0.01, 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.55, 0.6 , 0.65]

num_plots = len(ol_frac)



r_p = random.choice([0.5])  # pressure (GPa) of m/s segregation in reduced component    
o_p = random.choice([ 3])	# pressure (GPa) of m/s segregation in oxidised

r_c =  random.choice([1e-3])  # carbon content of reduced metals
o_c = random.choice([4])	  # carbon content of oxidised metals

r_s = random.choice([2]) 	# S content of reduced metals
o_s = random.choice([1])	# S content of oxidised metals





oxidised_fe = []
reduced_fe = []
oxfrac = []
finalW = []
finalNi = []
finalNi_conc =[]
finalNb =[]
finalCo = []
finalCo_conc =[]
finalCr = []
finalMo = []
finalV = []
finalSi = []



# bulk_earth is the CI_bulk Earth element content from the devolatilaised CI/pyrolite of, mainly, Sun &McDonough, 1995
# CI_total is CI composition from S&M , in ppm

if element == 'ni':
    CI_total = 1.05*10000
    CI_err = 100
    bulk_earth = 1.8*10000 
    
elif element == 'co':
    CI_total= 0.0500*10000
    CI_err = 50
    bulk_earth = 0.085*10000
    
elif element == 'w':
    CI_total= 0.0093*10000
    CI_err = 0.0093
    bulk_earth = 0.0171*10000 #12pppm in mantle from Konig 2011
    
elif element == 'nb':
    CI_total= 0.0240*10000
    CI_err = 30
    bulk_earth = 350
    
elif element == 'cr':
    CI_total= 0.26500*10000
    CI_err =500
    bulk_earth = 0.3883*10000
    
elif element == 'v':
    CI_total= 0.00560*10000
    CI_err =10
    bulk_earth = 96.6

elif element == 'mo':
    CI_total= 0.00090*10000 #ppb
    CI_err =50
    bulk_earth = 1.64
    
elif element == 'si':
    CI_total= 10.62*10000
    CI_err =1.62
    bulk_earth = 16
    
#elif element == 'mn':
#      CI_total= 0.1920*10000
#      CI_err =0.0192*10000
#      bulk_earth = 1760
    





     
    

df = pd.read_csv('merged_file.csv')

    


P = df['Pressure']
ol = df['silicate_olivine%']
c = df['C_metal']

fe = (df['fe_mant'])   
s = df['S_metal']   # S conc of metals

name= 'ave_'+element
name_err = element+'_err'



ol_all = df['silicate_olivine%'].unique()

df.loc[(df.C_metal >5), 'C_metal'] = 999

c_all = df['C_metal'].unique()
s_all = df['S_metal'].unique()
p_all = df['Pressure'].unique()

feall = fe.unique()

fe_ox_all =[12. ,
       14. , 16. , 18. , 20. , 22. , 24. ,  28. , 30. ]  #df['fe_mant'].unique()
# fe_red_all =  0.3 # use list [] if more than one - cant do twin axes tho
plot_number = 0

n_ox = 0


fig = plt.figure(figsize=(6, 4*num_plots))  # Adjust the figure size as needed

# Create a GridSpec with num_plots rows and 1 column
gs = gridspec.GridSpec(num_plots, 1)



for o_ol in ol_frac:
    oxidised_fe = []
    reduced_fe = []
    oxfrac = []
    finalW = []
    finalNi = []
    finalNb =[]
    finalCo = []
    finalCr = []
    finalMo = []
    finalMn = []
    finalV = []
    finalSi = []
    d_data =[]
    mean = []
    
    fin_el_ci = []; oxidised_d = []; reduced_d = []
    y_smooth_hi = 0; y_smooth_lo = 0 ; y_smooth_mean = 0
    out_el = 0
    
    
    for i in range (1,100):
        
        r_fe =   0.3   #  random.choice(fe_red_all)   amount of silicate Fe of the reduced body
        
        
        
        # select olivine amounts reduced, oxidised set at line 85
        # reduced olivine fraction set here
        r_ol = 0.3
       
        
        r_row  = (df.loc[(c==r_c) & (ol == r_ol)  & ( fe == r_fe) 
                    & (P==r_p) & (s == r_s)])
        
        r_el_raw = float(r_row[name].iloc[0])
        r_el_err_raw = float(r_row[name_err].iloc[0])
        
        """
        over all errors dominated by  the reduced metal silicate partitioning errors
        which can be very big for elements that are highly siderophile in reduced component 
        hence, restrict the overall error to 10% of mean if required 
        
        
        in both oxidised and reduced, 
        D cannot be less that zero - if it is when sampling, make it megligeable
        """
        if r_el_err_raw > (r_el_raw/10):
            r_el_err_raw=r_el_raw/10
            
        
        
        
        
        r_el = np.random.normal(r_el_raw, r_el_err_raw)
        if r_el <=0:
            print ('hello', r_el, r_el_raw, r_el_err_raw)
            r_el = 1e-8
           
        
        
        
               
        r_cf = float(r_row['cf'].iloc[0])
        r_mf = 1-r_cf
        
        o_fe =  random.choice(fe_ox_all)
 
      
       
        
        o_row  = (df.loc[(c==o_c) & (ol == o_ol)  & (fe ==o_fe) 
                    & (P==o_p) & (s == o_s)])
        
        o_el_raw = float(o_row[name].iloc[0])
        o_el_err_raw = float(o_row[name_err].iloc[0])
        
        o_el = np.random.normal(o_el_raw, o_el_err_raw)
        
        
        
        if o_el <=0.:
           
            n_ox += 1
            o_el = 1e-8  # o_el_raw/100
       
        
# """
# o_cf is oxidised body's core fraction
# o_mf is oxidised body's  mantle fraction

# ox_frac is the amount of oxidised body need to be added to a reduced amount of given Fe content
# such that the final Earth contains the requisite 6.25% mantle Fe

# """        
        o_cf = float(o_row['cf'].iloc[0])
        o_mf = 1-o_cf
        
        ox_frac = (6.25-r_fe)/(o_fe-r_fe)
        
        final_mant_size = ox_frac*o_mf + r_mf*(1-ox_frac)
        final_core_size = r_cf*(1-ox_frac) + o_cf* ox_frac
        
        ox_sil_el = 1./(o_mf+o_cf*o_el)
        ox_met_el = 1./(o_mf/o_el + o_cf)
        
        red_sil_el = 1./(r_mf+r_cf*r_el)
        red_met_el = 1./(r_mf/r_el + r_cf)
        
        
        
        final_mant_el = (ox_frac*o_mf*ox_sil_el + (1-ox_frac)*r_mf*red_sil_el)/final_mant_size
        final_core_el = (ox_frac*o_cf*ox_met_el + (1-ox_frac)*r_cf*red_met_el)/final_core_size
        
        CI_bulk = np.random.normal(CI_total, CI_err)   
        
        final_d_el = final_core_el/final_mant_el
        
        
        if final_d_el < 0:
            print ('bugger')
            
            break
        oxidised_d.append(o_el)
        reduced_d.append(r_el)
        oxfrac.append(ox_frac)
        oxidised_fe.append(o_fe)
        
        CI_bulk = np.random.normal(CI_total, CI_err) 
        
        if element == 'ni':
              
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
           
            finalNi.append(d_conc )
            fin_el_ci.append(final_d_el)
            out_el= finalNi
            
        elif element == 'co':
            
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            #print('mant', final_mant_el*CI_bulk*(1.-lv), final_mant_el_conc, lv*CI_bulk)
            
            fin_el_ci.append(final_d_el)
            fin_el_ci.append(final_d_el)
            finalCo.append(d_conc)

            out_el = finalCo
            
        elif element == 'w':

            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            fin_el_ci.append(final_d_el)
            finalW.append(d_conc )
            out_el= finalW
            
            
        elif element == 'v':
                
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            fin_el_ci.append(final_d_el)
            finalV.append(d_conc )
            out_el = finalV
        
        elif element == 'cr':
           
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            
            fin_el_ci.append(final_d_el)
            finalCr.append(d_conc )
            out_el = finalCr
            
            
        elif element == 'mo':
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            
            fin_el_ci.append(final_d_el)
            finalMo.append(d_conc )
            out_el = finalMo
            
        elif element == 'mn':
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            
            fin_el_ci.append(final_d_el)
            finalMn.append(d_conc )
            out_el = finalMn
        
        elif element == 'nb':
              
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            fin_el_ci.append(final_d_el)
            finalNb.append(d_conc )
            out_el= finalNb
          
        elif element == 'si':
            final_mant_el_conc = final_mant_el*CI_bulk*(1.-lv)+lv*CI_bulk
            final_core_el_conc = final_core_el*CI_bulk
            d_conc = final_core_el_conc/final_mant_el_conc
            
            fin_el_ci.append(final_d_el)
            finalSi.append(d_conc )
            out_el = finalSi
    
    
    d_data = pd.DataFrame(list(zip(oxfrac, oxidised_fe, out_el, fin_el_ci, reduced_d,oxidised_d, )),
                      columns = ['ox_frac', 'oxfe', element, 'd_before_LV', 'reduced_d', 'oxidised_d' ])
    
    #set y-limits for plots 
    if element == 'ni' or element =='co':
        ymin, ymax = 1, 2.5
    elif element == 'nb' or element == 'cr' or element =='v' or element =='mn':
        ymin, ymax = -3, 3
    elif element == 'w'  :
        ymin, ymax = 0, 3.0
    elif element == 'mo':
        ymin, ymax = 0, 3
    elif element == 'si':
        ymin, ymax = -3, 1
   
    
    x = d_data['oxfe']
    y = np.log10(d_data[element])
    #set x-limits for plots
    xmin, xmax = x.min(),x.max() 
    
    ex_x = xmin, xmax
    ex_y= ymin, ymax
    
    
   # origin = 'lower'
    extent = (ex_x + ex_y)
    
    mean = d_data.groupby(['oxfe']).mean()
    mean.reset_index(inplace=True)
    
    stdev = d_data.groupby(['oxfe']).std()
    stdev.reset_index(inplace=True)
    
    x_low = np.min(mean['oxfe'])
    x_hi = np.max(mean['oxfe'])
    x_new = np.linspace(x_low, x_hi,500)
    
    if element == 'nb':
        sig = 2.
    else:
        sig = 2.
             
    y_m = interp1d(mean['oxfe'],(mean[element]), kind = 'slinear', fill_value = 'extrapolate')   
    y_up = interp1d(mean['oxfe'], (mean[element]+sig*stdev[element]), kind = 'slinear', fill_value = 'extrapolate')
    y_lo = interp1d(mean['oxfe'], (mean[element]-sig*stdev[element]), kind = 'slinear', fill_value = 'extrapolate')
     
    y_smooth_mean=y_m(x_new)
    y_smooth_hi=y_up(x_new)
    y_smooth_lo=y_lo(x_new)
    
    y_smooth_lo[y_smooth_lo < 0] =1e-6
    
    
   
    ax1 = fig.add_subplot(gs[plot_number])
    
    if o_c == 1e-3:
        o_c = 0
    if o_s == 1e-3:
        o_s = 0
        
    if r_c == 1e-3:
        r_c = 0
    if r_s == 1e-3:
        r_s = 0
    
    if element == 'ni':
        
           # ax1.imshow(np.rot90(f), cmap='Blues', aspect = 'auto', extent = extent,  origin = 'upper', alpha = 0.5)   
        
            rect_ni = patches.Rectangle((9.,1.39), 30, 0.1, linewidth = 1, edgecolor = 'none', facecolor = 'darkred',
                                        alpha = 0.3)
            
                                
            ax1.add_patch(rect_ni)
        
            ax1.plot(x_new, np.log10(y_smooth_mean), 'b-.')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), alpha=.3, linewidth=.2, edgecolor = 'r')
           
            #  the following fills between the lines with a gradient but not centered around the mean, so its obvs crap. obvs. Looked nice tho....
            
            # polygon = ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), color = 'none', linewidth=0)
            # verts = np.vstack([p.vertices for p in polygon.get_paths()])
            # ymin, ymax  = verts[:, 1].min(), verts[:, 1].max()
            
            # gradient = ax1.imshow(np.array([np.interp(np.linspace(ymin,ymax, 200), [y1i, y2i], np.arange(2))
            #                     for y1i, y2i  in zip(np.log10(y_smooth_hi), np.log10(y_smooth_lo))]).T,
            #           cmap='turbo', aspect='auto', origin='lower', extent=[x.min(), x.max(), ymin, ymax])
            
            
            
            # gradient.set_clip_path(polygon.get_paths()[0], transform=plt.gca().transData)
            
            # ax1.plot(d_data['oxfe'], np.log10(d_data[element]), 'r.', alpha=0.4)
            #ax1.imshow(np.rot90(f), cmap='Blues', aspect = 'auto', extent = extent,  origin = 'upper')
          
            ax1.set_ylabel(r'Log D$_{10}$ Ni $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate('Ni', xy = (12.25,1.5), xytext=(13.,2.25), size = 24)
            ax1.annotate(str(o_ol*100) + '% Olivine$_{ox}$', xy = (12.25,2.2), xytext=(12.25,1.85), size = 10)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 1.19),  xytext=(12.25,1.25), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 1.15),  xytext=(12.25,1.15), size = 10)
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,1.2), xytext=(12.25,1.65), size = 10)

            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,1.05),  xytext=(12.25,1.04), size = 11)
            ax1.set_ylim([ymin,ymax])
            levels = [ 0.4, 0.6, 0.9]
    elif element =='co':
            
            ax1.plot(x_new, np.log10(y_smooth_mean), 'r-.')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'salmon', alpha=.3, linewidth=.5, edgecolor = 'b')
            
            rect_co = patches.Rectangle((0.,1.36), 30., 0.05, linewidth = 1, edgecolor = 'none', facecolor = 'navy',
                                            alpha = 0.3)
            ax1.add_patch(rect_co)   
            ax1.set_ylabel(r'Log D$_{10}$ Ni $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate('Co', xy = (12.25,1.5), xytext=(13.,2.25), size = 24)
            ax1.annotate(str(o_ol*100) + '% Olivine', xy = (12.25,2.2), xytext=(12.25,2.07), size = 10)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 1.9),  xytext=(12.25,1.97), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 1.8),  xytext=(12.25,1.87), size = 10)
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,1.2), xytext=(12.25, 1.17), size = 10)

            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,1.05),  xytext=(12.25,1.05), size = 11)
            ax1.set_ylim([ymin,ymax])
            levels = [ 0.4, 0.6, 0.9]
            
       
    elif element =='w':
            # from wang et al: The elemental abundances (with uncertainties) of the most Earth-like planet  
          
            rect_w = patches.Rectangle((0,1.38), 30., 0.39, linewidth = 1, edgecolor = 'none', facecolor = 'orangered',
                                            alpha = 0.40)
           
            ax1.add_patch(rect_w) 
            
            ax1.plot(x_new, np.log10(y_smooth_mean), color = 'orangered', ls = '-.')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), alpha=.3, linewidth=.2, edgecolor = 'r')
                  
            ax1.set_ylabel(r'Log D$_{10}$ W $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate('W', xy = (12.25,2.5), xytext=(13.,2.5), size = 24)
            ax1.annotate(str(o_ol*100) + '% Olivine ', xy = (12.25,1.05), xytext=(12.25,1.1), size = 11)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 0.7),  xytext=(12.25,0.92), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 0.7),  xytext=(12.25,0.73), size = 10)

            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,0.7), xytext=(12.25,0.35), size = 11)

            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,1.05),  xytext=(12.25,0.17), size = 11)
            ax1.set_ylim([ymin,ymax])
            levels = [ 0.4, 0.6, 0.9]
            
    elif element =='mo':
           
            rect_mo = patches.Rectangle((0,1.81), 30, 0.43, linewidth = 1, edgecolor = 'none', facecolor = 'darkred',
                                             alpha = 0.2)
            ax1.add_patch(rect_mo)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), color = 'orangered', ls = '-.')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), alpha=.3, linewidth=.2, edgecolor = 'r')
           
            
          
            ax1.set_ylabel(r'Log D$_{10}$ W $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate('Mo', xy = (12.25,2.5), xytext=(13.,2.5), size = 24)
            ax1.annotate(str(o_ol*100) + '% Olivine ', xy = (12.25,1.05), xytext=(12.25,1.1), size = 11)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 0.7),  xytext=(12.25,0.92), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 0.7),  xytext=(12.25,0.73), size = 10)

            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,0.7), xytext=(12.25,0.35), size = 11)

            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,1.05),  xytext=(12.25,0.17), size = 11)
            ax1.set_ylim([ymin,ymax])
        
    elif element == 'v':
            # from wang et al: The elemental abundances (with uncertainties) of the most Earth-like planet 
          
            rect_v = patches.Rectangle((0.,-0.05), 30, 0.3, linewidth = 1, edgecolor = 'none', facecolor = 'deepskyblue',
                                            alpha = 0.4)
            ax1.add_patch(rect_v)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), '-.', color = 'goldenrod')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'wheat', alpha=.5, linewidth=.2, edgecolor = 'goldenrod')
            
           
            ax1.annotate('V', xy = (12.25,0.8), xytext=(28.,0.75), size = 24)
            
            ax1.set_ylabel(r'Log D$_{10}$ V $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate(str(o_ol*100) + '% Olivine, ', xy = (12.25,0.0), xytext=(12.25,-0.26), size = 10)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.25,-0.38), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.25,-0.51), size = 10)
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,0.75), xytext=(12.25,0.7), size = 10)
            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,0.13),  xytext=(12.25,-.63), size = 10)
           
            ax1.set_ylim([-0.75,1.])
            levels = [ 0.4, 0.6, 0.9]
            
    elif element == 'mn':
        
            rect_mn = patches.Rectangle((0.,-0.47), 30, 1.25, linewidth = 1, edgecolor = 'none', facecolor = 'deepskyblue',
                                            alpha = 0.4)
            ax1.add_patch(rect_mn)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), '-.', color = 'goldenrod')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'wheat', alpha=.3, linewidth=.2, edgecolor = 'goldenrod')
            ax1.annotate('Mn', xy = (12.25,0.8), xytext=(13.,0.75), size = 24)
            
            ax1.set_ylabel(r'Log D$_{10}$ Mn $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate(str(o_ol*100) + '% Olivine, '+str(o_c)+' wt% C$_{core}$', xy = (12.25,0.0), xytext=(12.45,-0.26), size = 11)
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,-0.04), xytext=(12.45,-0.48), size = 11)
            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,0.13),  xytext=(12.45,-.6), size = 11)
           
            ax1.set_ylim([-0.75,1.])
            levels = [ 0.4, 0.6, 0.9]
            
    elif element =='cr':
            # ax1.imshow(np.rot90(f), cmap='BuGn', aspect = 'auto', extent = extent,  origin = 'upper')
            
            rect_cr = patches.Rectangle((0,0.39), 30, 0.15, linewidth = 1, edgecolor = 'none', facecolor = 'navy',
                                            alpha = 0.25)
            ax1.add_patch(rect_cr)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), '-.', color = 'darkmagenta')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'salmon', alpha=.3, linewidth=.2, edgecolor = 'darkmagenta')
            
            ax1.annotate('Cr', xy = (12.25,0.8), xytext=(27.,0.8), size = 24)
            
            ax1.set_ylabel(r'Log D$_{10}$ Cr $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate(str(o_ol*100) + '% Olivine,', xy = (12.25,0.2), xytext=(12.25,0.15), size = 10)
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.25,0.07), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.25,-0.0), size = 10)
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,0.01), xytext=(12.25,0.74), size = 10)
            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,0.13),  xytext=(12.25,-.09), size = 11)
            
            
            ax1.set_ylim([-0.2,1.])
            levels =  [ 0.4, 0.6, 0.9]
            
    elif element =='nb':
    # # Nb labelling
           # ax1.imshow(np.rot90(f), cmap='PuBuGn', aspect = 'auto', extent = extent,  origin = 'upper')
            ax1.axhline(y=-.1, linestyle = 'dashed', linewidth = 3, color='navy', alpha=0.7)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), '-.', color = 'teal')
            # assume the error is +/- 0.5 log unit given disparity in published regressions
            # ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'aquamarine', alpha=.3, linewidth=.2, edgecolor = 'g')
            ax1.fill_between(x_new,np.log10(y_smooth_mean)+0.5, np.log10(y_smooth_mean)-0.5, fc = 'aquamarine', alpha=.3, linewidth=.2, edgecolor = 'g')
            
            ax1.set_ylabel(r'Log D$_{10}$ Nb $^{[core]}/_{[mantle]}$', fontdict = font, fontsize = 14)
            ax1.annotate('Nb', xy = (12.25,2.15), xytext=(28,2.15), size = 24)
            ax1.annotate(r'Upper $D$ Limit', xy = (12.25,2.15), xytext=(25.2,-0.5), size =11)
            ax1.annotate(str(o_ol*100) + '% Olivine,', xy = (12.45,0.2), xytext=(12.45,-1.3), size = 10)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.45,-1.8), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, 0),  xytext=(12.45,-2.25), size = 10)
           
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.45,0.01), xytext=(12.45,1.8), size = 10)
            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,0.13),  xytext=(12.45,-2.7), size = 10)
            
            
            ax1.set_ylim([ymin,ymax])
           
            
    elif element =='si':
            # ax1.imshow(np.rot90(f), cmap='Oranges', aspect = 'auto', extent = extent,  origin = 'upper')
            rect_si = patches.Rectangle((0,-0.84), 30, 0.42, linewidth = 1, edgecolor = 'none', facecolor = 'navy',
                                            alpha = 0.3)
            ax1.add_patch(rect_si)
            
            ax1.plot(x_new, np.log10(y_smooth_mean), '-.', color = 'teal')
            ax1.fill_between(x_new,np.log10(y_smooth_hi), np.log10(y_smooth_lo), fc = 'aquamarine', alpha=.7, linewidth=.2, edgecolor = 'g')
            
            ax1.annotate('Si', xy = (12.25,0.6), xytext=(13.,0.3), size = 24)
            
            ax1.annotate(str(o_ol*100) + '% Olivine$_{ox}$', xy = (12.25,2.2), xytext=(12.25,1.85), size = 10)
            
            ax1.annotate('Red '+ str(r_s)+'wt% S$_{core}$; '+str(r_c)+' wt% C$_{core}$', xy = (12.25, -1.19),  xytext=(12.25,-1.65), size = 10)
            ax1.annotate('Ox: '+ str(o_s)+'wt% S$_{core}$; '+str(o_c)+' wt% C$_{core}$', xy = (12.25, -1.15),  xytext=(12.25,-1.9), size = 10)
            
            ax1.annotate(str(o_p) + ' GPa  - Peak Pressure in'+'\n'+'             oxidised material', xy = (12.25,-2.2), xytext=(12.25,-2.45), size = 10)

            ax1.annotate(str(lv*100) + ' % late veneer', xy = (12.25,0.13),  xytext=(12.45,-2.8), size = 11)
            
            ax1.set_ylim([-3,1.])
            levels =  [ 0.4, 0.6, 0.9]
    
    
    
    ax1.set_xlabel(r'Fe$_{(in}$ $_{oxidised}$ $_{mantle)}$ wt%',fontdict=font)
    ax1.set_xlim(ex_x)
    
    
    ax1.tick_params(axis='y', left=True, labelleft=True)
    
    #  mark chondrite/meteorite groups
    

    
    ax2 = ax1.twiny()
    
    #  sort out top axis - this is a faff but its here for future reference
    #  actually, thats the only reason. Missed match of the day doing this.... Kind of wish I hadn't bothered...

    
    new_ticks_list = []
    inv=[]
    ox_fe_ticks = [.5, .4, 0.3, 0.25,  .2]
    
    for o in ox_fe_ticks:
               
          oxf = (6.25-(1.-o)*r_fe)/o
          
          checksum = oxf*o+r_fe*(1-o)
          print (r_fe, o, oxf, checksum)
          new_ticks_list.append(oxf)
             
            
    new_ticks = np.array(new_ticks_list)  
  
    for tick_x in new_ticks:
        inv.append(1/tick_x) 
    
    def func(fx, a, c):
        y = a*1/fx + c
        return y
    
    popt, pcov = curve_fit(func,  new_ticks, ox_fe_ticks)
    
    def tick_function(X, popt):
        V = popt[0]/(X-popt[1])
        return ["%.1f" % z for z in V]
    
    
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_ticks)
    ax2.set_xticklabels(tick_function(new_ticks, popt))
    ax2.set_xlabel(r'Fraction of oxidised body', fontdict=font, fontsize=12)
     
    new_y_ticks_list = []
     
    
    if element =='ni':
        conc_ticks = [4000, 1900, 1000, 500]
        
    elif element =='co':
        conc_ticks = [ 140,100, 75, 50, 30, 20]
        
    elif element =='w':
        conc_ticks = [ 200, 50, 20, 5, ]
        
    elif element =='nb':
        conc_ticks = [5, 50, 100, 300, 400, 500]
    
    elif element =='cr':
        conc_ticks = [ 3500, 2500, 2000, 1500]
    
    elif element =='v':
        conc_ticks = [  120,85,  50, 30 ]
    
    elif element =='mo':
        conc_ticks = [  .5, 0.05]
    
    elif element =='si':
        conc_ticks = [ 10, 18, 21, 22, 23]
        
    #elif element =='mn':
    #    conc_ticks = [  1000, 1500, 2000, 2300]
        
        
    cf =0.315    
    for e in conc_ticks:
        
        equivalent_D = np.log10((bulk_earth/(e*cf))-1/cf+1)                        
        new_y_ticks_list.append(equivalent_D)
        new_y_ticks = np.array(new_y_ticks_list)  
    
        
    #
    print ('new y ticks', new_y_ticks)
    
    
    ax3 = ax1.twinx()
    ax3.set_ylim(ax1.get_ylim())
    ax3.set_yticks(new_y_ticks)
    ax3.set_yticklabels(conc_ticks)
    
    
    if element =='w' or element =='mo':
        ax3.set_ylabel (r'Mantle concentration (ppb)', fontdict=font, fontsize = 11)
        
        #  mark chondrite/meteorite groups
        r_chond = patches.Rectangle((24.7,2.6), 0.5, 2, linewidth = 1, edgecolor = 'none', facecolor = 'darkviolet',
                                    alpha = 0.3)
        ax1.add_patch(r_chond)
        ax1.annotate(r'$R$', xy = (24.7,3.0), xytext=(24.7,2.65), size = 12)
    elif element =='si':
        ax3.set_ylabel (r'Mantle concentration (wt%)', fontdict=font, fontsize = 11)
        r_chond = patches.Rectangle((24.7,2.6), 0.5, 2, linewidth = 1, edgecolor = 'none', facecolor = 'darkviolet',
                                    alpha = 0.3)
        ax1.add_patch(r_chond)
        ax1.annotate(r'$R$', xy = (24.7,3.0), xytext=(24.7,2.65), size = 12)
    elif element =='cr' or element =='v':
        ax3.set_ylabel (r'Mantle concentration (wt%)', fontdict=font, fontsize = 11)
        r_chond = patches.Rectangle((24.7,0.85), 0.5, 0.15, linewidth = 1, edgecolor = 'none', facecolor = 'darkviolet',
                                    alpha = 0.3)
        ax1.add_patch(r_chond)
        ax1.annotate(r'$R$', xy = (24.7,1.0), xytext=(24.7,0.875), size = 12)
    else:
        ax3.set_ylabel (r'Mantle concentration (ppm)', fontdict=font, fontsize = 11)
        
        #  mark chondrite/meteorite groups
        r_chond = patches.Rectangle((24.7,2.3), 0.5, 2, linewidth = 1, edgecolor = 'none', facecolor = 'darkviolet',
                                    alpha = 0.3)
        ax1.add_patch(r_chond)
        ax1.annotate(r'$R$', xy = (24.7,2.3), xytext=(24.7,2.32), size = 12)
    
    plot_number = plot_number + 1
    
plt.tight_layout()
# save_results_to = '/Users/jonw/Google Drive/My Drive/Low_P_core_model/monte_plots_output/'
# plt.savefig(save_results_to+element+ '_'+ str(o_c)+'_%oxC_'+ str(r_c)+'_%redC_'+str(lv*100)+'_lateV_' + str(o_s)+'_%oS_'+ str(r_s)+'%rS_'+str(lv*100)+'%l_v_sept.pdf', dpi = 300)
#

plt.show()
