# accrete-o-matit
a script for tracking mantle/core partitioning 


This is a script for tracking the core mantle partitioning behaviour of elements during core formation.
It allows for core/mantle partitioning to occur in the presence of an olivine fraction in the silicate
and variable C and S contents in the metal (these are set by the script.
The script uses a the metal/silicate experimental data from the literature (see paper for deets) 
and olivine melt partitioning data, with the metal activities using the epsilon formula of Ma 
and interaction parameters from a variety of sources.  Full refs are given in the paper and I'll add here when 
I get a chance.  

A warning - script is a bit messy  and comments within may not be useful. Theres also bits in the script that are there for other uses.  Or because I thought they were a good idea at the time....

another warning - its not exactly efficient!  because it loops over the various parameters (P, Mantle Fe contnets,  metal C & S contents and melt Olivine contenst) and uses the python 'uncertainties' package and is sloooooooww. To make it even vaguely acceptable I split the different peak pressures of accretion (here 0.5 and 3 GPa)  and run two (or more) instances. If this isn't done, it would be quicker to accrete your own planet and measure the partitioning directly.  as it is, takes about 14 hours to do an oxidised run on my slighlty battered M1 MacBook. 

In the folder labelled 'Output and plotting' theres a script that will merge two csv files for use in the script that plots out the graphs. 

I will *obviously* edit this readme so it makes sense to someone other than me.  
