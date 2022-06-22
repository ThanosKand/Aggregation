# Aggregation

Model: The basic model is run by the call

simo = coagfun(a,alpha,epsilon,Ptotal,[],Tmax,seasonal)

where the input parameters are
a	the self-similarity parameter e.g. 1.7 to 2.1
alpha	the stickiness factor, in the range 0 to 1.
epsilon	turbulent dissipation rate e.g. 1E-4 units: m2 s-3
Ptotal	rate of production of primary particles e.g. 1E5 units [µg C m-2 day-1]
[]	is a space holder for an input structure – not implemented here
Tmax	the time span of simulations e.g. 3*365 units [days]
seasonal	0 for steady and 1 for a simple seasonal production cycle.

All relevant output is passed through the structure simo. Most relevant are
simo.N	Aggregate distribution at Tmax:  Size and density resolved number [m-3] 
simo.M	Aggregate distribution at Tmax:  Size and density resolved mass [µg C m-3]  
simo.p	All run parameters
simo.w	Sinking speed [m day-1]
simo.d	Excess density [µg µm-3]
simo.m	Mass [µg]

For the sensitivity plots, sweeps across parameter space were set up using the batchrun script containing instruction sets such as

clear all

sima17 = coagfun(1.7,0.5,1E-4,1E5,[],1*365,0);
sima18 = coagfun(1.8,0.5,1E-4,1E5,[],1*365,0);
sima19 = coagfun(1.9,0.5,1E-4,1E5,[],1*365,0);
sima20 = coagfun(2.0,0.5,1E-4,1E5,[],1*365,0);
sima21 = coagfun(2.1,0.5,1E-4,1E5,[],1*365,0); 
sima22 = coagfun(2.2,0.5,1E-4,1E5,[],1*365,0);
sima23 = coagfun(2.3,0.5,1E-4,1E5,[],1*365,0);

save sima1y0;

This performs 7 model runs, each for a stickiness parameter of 0.5, turbulent dissipation rate of 10-4 m3 s-3, a total productivity of 1E5 µg C m-2 day-1 and over a 1 year time span with constant (i.e. non-seasonal) forcing. The only difference being the value of the self-similarity parameter; ranging here from 1.7 to 2.3. All out put structures are save in the file sima1y0.mat that can be reloaded for subsequent analysis.

To plot sensitivity results, use batchplotx. Reset the paramset to ‘a’, ‘e’ or ‘s’ for the appropriate sensitivity sweeps. These produce the panels found in Figure 3 of the manuscript. There are additional panels that present mean excess density and sinking speed within size bins of the emerging aggregate community, as well as an indication of the sinking speed at the maximum flux. 

