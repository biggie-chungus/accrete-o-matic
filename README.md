# accrete-o-matic
a script for tracking mantle/core partitioning 


**accretion_1.py**
This is a script for tracking the core mantle partitioning behaviour of elements during core formation.

It allows for core/mantle partitioning to occur in the presence of an olivine fraction in the silicateand variable C and S contents in the metal (these are set by the script. Following on from other work, the script increments the pressure of core/mantle equilibrium to a final value (peak pressure - variable 'p_splat' on line 720).  

The script uses a the metal/silicate experimental data from the literature (see paper for deets) and olivine melt partitioning data, with the metal activities using the epsilon formula of Ma and interaction parameters from a variety of sources.  Full refs are given in the paper and I'll add here when I get a chance.  The partitioning data for Ni and Co should only be used 0 to 5 GPa (range over which it were fit).  Equations can obviously be changed and pressures likewise.  Note: some metal/silicate regressions in the literature use elements concentrations rather than mole fractions.  Additionally, 1/T term for W has been modified to include an offset from the Jennings et al paper which had smaller value at NBO/T's of ~2.7. Given there is a large uncertainty in the compositional effects of silicate melt on the M/S partitioning behaviour of W, which should be especially profound in the Fe-rich silicate melts here.

A warning - script is a bit messy  and comments within may not be useful. Theres also bits in the script that are there for other uses.  Or because I thought they were a good idea at the time....

another warning - its not exactly efficient!  because it loops over the various parameters (P, Mantle Fe contnets,  metal C & S contents and melt Olivine contenst) and uses the python 'uncertainties' package and is sloooooooww. To make it even vaguely acceptable I split the different peak pressures of accretion (here 0.5 and 3 GPa)  and run two (or more) instances. If this isn't done, it would be quicker to accrete your own planet and measure the partitioning directly.  as it is, takes about 14 hours to do an oxidised run on my slighlty battered M1 MacBook. 
**remember to make sure each instance outputs to a file with a different file name (line 688)**


**Output_and_plotting_scripts **folder

When the accretion_1.py script has done its thing, put the output into this folder and run the 'merge_files' script.
This will join the output for different runs, if required, with the merged_file used by the plotting script.

In the folder labelled 'Output and plotting' theres a script that will merge two csv files for use in the script that plots out the graphs. 

I will *obviously* edit this readme so it makes sense to someone other than me and even then.....

Jon Wade. 
