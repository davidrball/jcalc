# jcalc

Code to read in data from runs from the HARM simulation, perform coordinate transformations and appropriate conversions from 
code units to physical units.  Also includes a function for covariantly calculating the current density in the Kerr metric
given a magnetic field profile.

things to include:
displacement current (system is fairly relativistic so this term likely matters),
put the del(E) term in the proper units,
output coordinate invariant, j_nu j^nu rather than the individual components.
