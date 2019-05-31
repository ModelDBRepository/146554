Version 1.1 November 2nd, 2012

--------------------------------------------------------------------------

Hidden state and parameter estimation for a network model of sleep:

The included set of files reproduce figures in:

Sedigh-Sarvestani, Schiff, Gluckman,'Reconstructing Mammalian Sleep
Dynamics with Data Assimiliation', PLoS Comp Biol, 2012 (In
press). DOI:10.1371/journal.pcbi.1002788

The data assimilation framework uses the Diniz Behn (DB) and the
Fleshner, Booth, Forger, Diniz Behn (FBFD) models of sleep: DB model
from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.  FBFD
model from:Fleshner, Booth, Forger, Diniz Behn, Philos Transact A Math
Phys Eng Sci 369(1952):3855-83, 2011.

Implemented by Madineh Sedigh-Sarvestani, Steven Schiff and Bruce
Gluckman.  Contact Madineh for questions.
(m.sedigh.sarvestani@gmail.com).

--------------------------------------------------------------------------

The .zip file containes a folder for each figure and a folder to
generate data from the DB and FBFD model which serves as the true data
as well as the noise-added 'observations' throughout the figures.

In the folder for each figure is a standalone .m file which will
produce a figure similar to that in the paper. Many figure files share
identical sub-functions, such as the UKF function, but these have been
reproduced in each figure .m file for completeness.

Parameter values are stored outside of these folders to assure
universality across files.

CI_eps.mat contains CI_eps, the default covariance inflation
multiplier (gets multiplied by variance of variable) that gets applied
for all UKF steps which use default CI matrix.

Note: Figures 3 4 and 6 contain reconstruction of several days worth
of data at 0.5 dT seconds sampling rate and therefore take longer than
other figures to run, on the order of several hours.

--------------------------------------------------------------------------

Brief description of each figure (see Sedigh-Sarvestani et al. for
detailed methods and interpretation of results):

Figure1: generates data from DB and FBFD models of sleep and plots
short (DB) and long (FBFD) cycles.

The DB model of sleep contains 5 nodes, two wake_active, one
NREM-active, and one REM-active. Each node is described by its firing
rate and output neurotransmitter concentration.

The FBFD model of sleep is similar to the DB model, with the addition
of the SCN as an additional node. The SCN causes the output of this
model to have 24 hour periodicity entrained to the light-dark cycle.

Figure 2: Uses noise-added F-LC (firing rate of wake-active LC
variable) as an observable and reconstructs the remaining 11 variables
of the DB model. Panels A and B show that reconstruction can be
improved as a function of UKF covariance inflation.

Figure 3: Creates the empirical observability (EOC) matrix. EOC(x,y)
is a metric of how well variable x is reconstructed by observation of
variable y (via UKF). The EOC can be used to gain intuition regarding
the partial observability properties of the DB model. It can also be
used (as we show in Figure 4) to optimize covariance inflation values
to improve reconstruction.

Figure 4: Optimizes covariance inflation parameter for variables delta
and F_R (firing rate of REM-active group). Only these variables are
considered because optimization of CI for other variables does not
significantly improve reconstruction performance. Once optimized
values are determined, the EOC for default, optimal CI_delta, and
optimal CI_delat and F_R is plotted. Improvements in reconstruction
after optimization of CI can be seen by comparing these EOC matrices.

Figure 5: This is the first parameter estimation figure. The UKF model
contains an arbitrary value for the parameter gALC that is different
from the parameter used in the model that generated the data. UKF
reconstruction is carried out iteratively with parameter estimation in
30 minute long windows. The method of parameter estimation is similar
to a 'shooting method' approach. In each parameter estimation step, a
group of test-trajectories are projected forward where each
test-trajectry is generated from a model carrying a slightly different
parameter. The trajectory with the smallest distance to xhat
(reconstructed UKF estimate from last iteration) determines which
parameter gets passed on for the next UKF iteration. In order to
prevent divergence of test-trajectories, each trajectory is initilized
at xhat at 1 minute intervals.

Figure 6: This is the second parameter estimation figure. This is also
the first figure wherein the UKF model (modified DB) significantly
differs from the model that generated the observations (FBFD in this
case).  The DB model is slightly modified to include an additional
input parameter pSCN into the sleep-wake nodes. The FBFD model, from
which observations are generated, includes an SCN node with its own
firing rate and neurotransmitter that results in 24 hour periodicity
in the sleep-wake pattern. Using the same shooting-method parameter
esimation approach, pSCN is estimated. Thus, even though the UK filter
model lacked specific components of the underlying system (SCN
periodicity), it was nonetheless able to accurately estimate this
rhythm.

Figure 7: This figure is tied to figure 8. Figure 7 produces
observations by mapping overall sleep-state of the generated data set
to first and second order statistics of the firing rate variables
during each sleep-state. In this way, mapped firing rates for F_LC,
F_DR, F_VLPO, and F_R are generated and used as variables in the UKF
to reconstruct the remaining hidden variables. Although these
SOV_mapped observables clearly lack certain details as compared to
direct observation of firing rates, UKF reconstruction remains
reasonably accurate.

Figure 8: Figure 8 uses the observations generated in figure 7 to
estimate unknown parameter gALC, much in the same way as was done for
figure 5.
