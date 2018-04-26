
# coding: utf-8

# # Statistical constraints on observing low-angle normal fault seismicity

# ## Introduction
# 
# Low-angle normal faults (LANFs), with dips less than 30$^\circ$, are
# well described in the geologic record. They are hypothesized to have
# an important role in accommodating large-magnitude continental
# extension \citep{howard1987crustal} and crustal thinning
# \citep{lister1986detachment}, and their recognition has been one of
# the most important developments in tectonics over the past several
# decades \citep{wernicke2009detachment}. However, despite widespread
# field observations of inactive LANFs, and their central role in modern
# extensional tectonic theory, they remain enigmatic and contentious
# structures, and it is not clear if they are seismically active at low dip 
# angles in the upper crust. This is for two reasons: because brittle faulting
# on LANFs is in apparent conflict with standard Andersonian rock mechanical
# theory as typically applied to the upper crust
# \citep{axen2004lanfmech}, and because observations of active faulting
# on LANFs is sparse and at times ambiguous \citep{wernicke1995seis}. A
# considerable amount of research has been performed to address the
# former concern, reconciling LANF slip with rock mechanics \citep [e.g.,]
# []{axenbartley1997, collettini2011lanfmech}. The latter issue, the paucity of 
# definitive LANF earthquakes in focal mechanism catalogs, has inhibited
# hypothesis testing of LANF fault theory, and has also contributed to a mode
# of thought where the  absence of evidence for LANF seismicity is taken
# as evidence that LANFs are inactive or aseismic \citep{jackson1987,
# collettinisibson2001}. However, the lack of observed seismic slip on
# continental LANFs may be simply due to the rarity of earthquakes on
# LANFs, coupled with the small number of potential active structures.
# 
# In this work, we choose to directly address the question of whether
# the lack of observed large earthquakes on LANFs may be 
# interpreted as an indication that
# LANFs do not slip seismically, or if it is better explained as an
# effect of the small number of active seismogenic LANFs with low moment
# accumulation rates (and hence long recurrence intervals). 
# We then calculate the maximum likelihood of observing a significant
# earthquake on a LANF by treating all potentially active LANFs described in the 
# literature as seismically active at their surface dip angles 
# throughout the upper crust. %, and displaying typical seismic behavior.
# Under these assumptions, we create synthetic earthquake catalogs with
# Gutenberg-Richter and `characteristic' frequency-magnitude distributions
# based on each fault's geometry and slip rate, and then we calculate the
# probability of observing earthquakes on each, and any, fault over different
# durations of observation.
# 

# ### LANF Slip, Mohr-Coulomb Failure Theory, and Seismicity

# Areas of the crust undergoing active extension are generally assumed
# to have a subvertical maximum compressive stress.  Mohr-Coulomb
# theory, as applied to the crust, predicts that a fault with a typical
# coefficient of friction for rocks (0.6--0.8) should lock up if it is
# oriented at an angle greater than 60$^\circ$ to the maximum compressive stress,
# and new, optimally oriented faults should form \citep{sibson1985}.  Therefore,
# for normal faults with dips less than 30$^\circ$, either much lower
# fault friction or elevated pore fluid pressure is required for fault slip.
# 
# Evidence for seismic slip on LANFs is sparse.  This is partly due to the
# ambiguity of the rupture plane in earthquake focal mechanisms, which have
# nodal planes 90$^\circ$ apart, and therefore have a
# high-angle nodal plane if they have a low-angle nodal plane.  Without
# ancillary information indicating which plane slipped, searches of earthquake
# catalogs cannot yield unique results as to whether they contain LANF events.
# Several collections of normal fault earthquakes
# with known surface breaks \citep{jackson1987, collettinisibson2001}, thereby
# resolving dip ambiguity, contain no low-angle events, though the collections
# are small ($\le$ 25 events).  Some candidate events exist, but they are
# undersea \citep[e.g.,][]{abers2001} or difficult to verify \citep[e.g.,]
# []{doser1987ancash}.

# ## Potentiall Active LANFs

# Over the past decade or so, many field studies have found evidence for LANF 
# activity in orogens throughout the world. These studies typically find arrays of 
# Quaternary normal fault scarps on the fault traces and/or in the hanging walls 
# of mapped or inferred low-angle detachment faults \citep [e.g.,]
# []{axen1999baja}. Some studies also have bedrock thermochronology data from the 
# exhumed footwalls of the detachments that is suggestive of ongoing rapid 
# exhumation \citep [e.g.,][]{sundell2013lunggar}, although this data does not 
# preclude a recent cessation of faulting. In some cases, additional evidence for 
# LANF activity comes from geophysical data such as GPS geodesy \citep [e.g.,]
# []{hreinsdottir2009altotib} and seismic waves \citep [e.g.,][]{doser1987ancash}.
# 
# We have compiled all potentially active LANFs with known subareal
# fault traces from a thorough review of the literature, finding twenty
# total (Figure~\ref{fig:lanf_map}).  We have then mapped the approximate fault
# traces into a GIS file (available at https://github.com/cossatot/LANF\_gis), 
# with metadata such as slip rate and source. Though the fault traces of many
# LANFs considered here are obscured by vegetation or agriculture, others
# display large fault scarps in Quaternary sediments, particularly those in the
# dry climates of Tibet \citep[e.g.,][]{styron2013slr, kapp2005nqtl} and the
# western US \citep[e.g.,][]{axen1999baja, hayman2003dv}, which are commonly
# interpreted as evidence for past seismic slip.  About half are in Tibet,
# consistent with hypotheses that LANFs and metamorphic core complexes
# form in areas of hot, thick crust \citep [e.g.,][]{buck1991mcc}.  The
# rest are distributed through other areas of active continental
# extension: the North American Basin and Range, the Malay Archipelago,
# western Turkey, Italy, and Peru. Several of the most-commonly cited
# candidates for seismically active LANFs were not included because they
# do not have a clearly-defined, mappable fault trace, which is
# necessary for our earthquake likelihood calculations.  These include the 
# submarine core complexes in the Woodlark Basin \citep{abers2001}, the fault
# responsible for the 1995 Aigion, Greece earthquake \citep{bernard1997}
# and other potential LANFs underneath the Gulf of Corinth, and the
# fault responsible for the 1952 Ancash, Peru earthquake
# \citep{doser1987ancash}.

# In[ ]:


# show map and figure caption


# ## Likelihood of observing a LANF event
# ### Earthquake likelihood on an individual LANF

# To estimate the likelihood of observing a significant earthquake on an
# individual LANF over some contiguous time window of length $t$ (in
# years), we perform a Monte Carlo simulation in which we create 4000
# synthetic time series of earthquakes, with unique values for fault
# geometry and slip rate for each time series. Then, for each time
# series we calculate the fraction of unique time windows of length $t$
# in which an earthquake as large or larger than a given magnitude
# occurs.  We take this value as the probability of observing an
# earthquake greater than or equal to moment magnitude $M$ over
# time period $t$, which we will refer to in general as $P(M,t)$.
# 
# The geometry for each fault is estimated based on the length of the
# fault trace, the dip of the fault, and the estimated seismogenic
# thickness or fault locking depth in the area.  The fault is treated as
# planar for simplicity of calculations, even though the exposed
# footwalls of many detachment faults are highly corrugated.  We
# determine the fault length by measuring the approximate length of the
# mapped fault trace perpendicular to the assumed extension direction;
# for faults that change dip significantly along strike (e.g., the Dixie
# Valley fault), we only consider the low-angle segments of the fault.
# Values for the dip are taken from the literature in most cases, and
# measured of the dip of footwall triangular facets (interpreted as the exhumed
# fault plane) from SRTM data otherwise. In all cases, ranges of fault
# geometries are considered, encompassing the degree to which the values are 
# known. The fault locking depth is assumed to be 10 km in the absence of other
# evidence (such as a geodetic study, \citep[e.g.,][]{hreinsdottir2009altotib}).
# 
# Slip rates of the 20 LANFs are similarly gathered from the
# literature if possible, or given broad ranges if not documented
# (e.g., 1--10 mm yr$^{-1}$).  In the Monte Carlo simulation, samples
# for slip rate and dip are drawn from uniform distributions defined
# by the maximum and minimum values.  Based on field observations,
# some faults have dip ranges that go above 30$^\circ$; this study
# only considers slip on faults shallower than this, so for these
# faults, dip values are sampled from the minimum to 30$^\circ$ and
# the resulting probabilities are then multiplied by the fraction of
# the dip range that is under 30$^\circ$.

# Each earthquake synthetic earthquake sequence is generated by
# randomly sampling either 50,000 events from a tapered Gutenberg-Richter (GR)
# distribution with corner magnitude $M_c = 7.64$ and $\beta = 0.65$
# (from values estimated by \citet{birdkagan2004f_m} for continental
# rifts), or a 25,000 events from `characteristic' distribution. It is not
# certain which distribution more appropriately describes seismicity on a single
# LANF, though studies of many individual fault rupture histories suggests
# that the characteristic distribution is more accurate \citep{hecker2013eqdist}.
# The smaller number of samples from the characteristic distribution is due to 
# the increased computation time associated with a higher proportion of large 
# events.  The samples are taken 
# from an interval $M = [5.0, \, M_{max}]$, where $M_{max}$ is calculated as the
# moment magnitude $M$ resulting from fault slip $D$ = 15 m over a fault of
# length $L$ cutting through a seismogenic thickness $z$ at dip $\delta$,
# given the relations
# 
# \begin{equation}
#  M_o = \mu L z D \,/ \, \sin \delta 
#  \end{equation}
# 
# and
# 
# \begin{equation}
# M = 2/3 \; \log_{10} (M_o) - 6
# \end{equation}
# 
# where $\mu = 30$ GPa is the shear modulus and $M_o$ is the seismic
# moment in N m \citep{kagan2003pepi}.  The characteristic distribution has a
# large-magnitude mode corresponding to $D$ = 1.5 m on the fault, a typical
# slip distance for normal fault events \citep[e.g.][]{wesnousky2008displacement}.
# The distributions are shown in Figure~\ref{fig:fms}.
# 
# These calculations rely on two important assumptions that warrant some
# discussion.  The first is that each earthquake ruptures the entire
# fault patch uniformly.  Though this is highly improbable fault behavior,
# the long-term statistical distribution of earthquake recurrence is 
# insensitive to assumptions about slip distribution in individual events:
# if $n$ different, equal fault patches rupture independently, each 
# requires $n$ times the interseismic strain accumulation time to rupture with
# some earthquake of magnitude $M$ versus a single fault rupturing uniformly; 
# but on the whole, magnitude $M$ events would happen with the same long-term
# frequency.  The next assumption is that earthquakes are ordered randomly and
# separated by the time necessary for sufficient strain to accumulate for each
# earthquake to occur.  This means that foreshock and aftershock sequences
# and other types of event clustering are not taken into account.  However,
# the modal inter-event times for earthquakes $\ge M \,6$ or so are greater than
# a hundred years for many LANFs [supplemental figure xx], so the ordering of
# events does not impact the results, as this is longer than our maximum
# observation window.  Furthermore, any clustering resulting in
# event spacing less than the observation window would decrease $P(M,t)$, and
# we are concerned with calculating the maximum $P(M,t)$.

# ### Setting up the problem

# #### Import necessary modules

# In[1]:


#%pylab inline
#%config InlineBackend.figure_format = 'retina'


# In[29]:

import matplotlib
matplotlib.use("Agg")
from matplotlib.pyplot import *


# In[3]:


import sys
sys.path.append('../eq_stats')

import numpy as np
import pandas as pd
import eq_stats as eqs
import time
from joblib import Parallel, delayed
from itertools import chain


# #### Read in fault data table
# 
# Makes a Pandas dataframe of fault data (length, slip rates, etc.)

# In[4]:


f = pd.read_csv('../data/lanf_stats.csv', index_col=0)


# #### Define some variables to be used later

# In[8]:


n_cores = -3
n_eq_samp = int(3e4) # number of earthquakes in time series
time_window = np.hstack( (1, np.arange(5, 105, step=5) ) ) # observation times
mc_iters = int(2e3) # number of Monte Carlo iterations
mc_index = np.arange(mc_iters, dtype='int')
mc_cols = ['dip', 'Ddot'] + [t for t in time_window]
max_eq_slip = 15 #m
char_eq_slip= 1.5 #m
Mc = 7.64 # Hypothetical corner magnitude for Continental Rift Boundaries


# Make list of minimum search magnitude $M_{min}$, and then make MultiIndex for Pandas dataframes

# In[9]:


min_M_list = [5, 5.5, 6, 6.5, 7, 7.5]

df_ind_tuples = [[i, M] for i in mc_index for M in min_M_list]
df_multi_ind = pd.MultiIndex.from_tuples(df_ind_tuples, names=['mc_iter','M'])

rec_int_bins = np.logspace(1, 5)


# ### Calculate $P(M,t)$ for Gutenberg-Richter frequency-magnitude distributions

# #### Define a function for Joblib Parallel to calculate probabilities for each iteration.
# 
# Function is defined here so it can access all variables generated by script, not just passed variables.  This makes the code cleaner even if it's not very abstracted.
# 
# Here is what this function does:
# 
# - Get the dip, Ddot and maximum earthquake magnitude for each iteration.
# - Take this info and make the earthquake sequence:
#     - Take the max earthquake magnitude and make a frequency-magnitude distribution based on a Gutenburg-Richter exponential model.
#     - Take 50k samples from this distribution, 
# - Make an earthquake time series form the EQ sequence
#     - Calculate the interseismic strain accumulation time for each event
#     - Separate each earthquake in the sequence with the appropriate number of years with no events.
# - Calculate the probability of observation
#     - Run a rolling maximum for each $t$ in [1, 5, 10, 15, ..., 95, 100]
#     - Calculate the observation probability above $M_{min}$ in [5, 5.5, 6, 6.5, 7, 7.5]
# - Calculate inter-event times for EQs $\ge \, M$

# In[7]:


def calc_iter_probs(iter):
    df_iter = fdf.loc[iter].copy()
    df_iter['dip'] = mc_d['dip_samp'][iter]
    df_iter['Ddot'] = mc_d['Ddot_samp'][iter]

    # Generate EQ sample/sequence from F(M) dist.
    m_vec = np.linspace(5, mc_d['max_M'][iter], num=1000)
    fm_vec = eqs.F(m_vec, Mc=Mc)
    M_samp = eqs.sample_from_pdf(m_vec, fm_vec, n_eq_samp)
    Mo_samp = eqs.calc_Mo_from_M(M_samp)
    
    # Make time series of earthquakes, including no eq years
    recur_int = eqs.calc_recurrence_interval(Mo=Mo_samp, 
                                             dip=mc_d['dip_samp'][iter],
                                             slip_rate=mc_d['Ddot_samp'][iter],
                                             L=params['L_km'],
                                             z=params['z_km'])

    cum_yrs = eqs.calc_cumulative_yrs(recur_int)
    eq_series = eqs.make_eq_time_series(M_samp, cum_yrs)
    
    # calculate probability of observing EQ in time_window
    for t in time_window:
        roll_max = pd.rolling_max(eq_series, t)
        df_iter[t] = (eqs.get_probability_above_value(roll_max, min_M_list)
                      * mc_d['dip_frac'])
        
    # calculate histgrams of recurrence intervals
    rec_int_counts_df = rec_int_df.loc[iter].copy()
    for mm in np.array(min_M_list):
        ints = np.diff( np.where(eq_series >= mm) )
        rec_int_counts_df.loc[mm] = np.histogram(ints, bins=rec_int_bins)[0]
    #print df_iter.as_matrix().shape
    return df_iter,rec_int_counts_df


# #### Iterate through the faults in the fault database, doing all the calculations for each.
# 
# The setup of this for loop is basically this:
# 
# - Make DataFrame for each fault.
#     - Columns are dip, Ddot, and observation time windows.
#     - Rows are values for each Monte Carlo iteration.  Values for time windows are calculated probabilities.
#     
# - Calculate maximum earthquake magnitude for each MC iteration.
# 
# - Run the above 'calc_iter_probs' function (parallelized over the MC iterations) and concatenate the results

# ### Note! Not actually running the calcs here, as it's a looong process!

# In[7]:

#for fault in list(f.index):
#    fdf = pd.DataFrame(index=df_multi_ind, columns=mc_cols, dtype='float')
#    params = f.loc[fault]
#    mc_d = {}
#    mc_d['dip_samp'], mc_d['dip_frac'] = eqs.dip_rand_samp( params['dip_deg'], 
#                                                         params['dip_err_deg'], 
#                                                         mc_iters)

#    mc_d['Ddot_samp'] = eqs.Ddot_rand_samp(params['slip_rate_mm_a'],
#                                           params['sr_err_mm_a'], mc_iters)
#
#    mc_d['max_Mo'] = eqs.calc_Mo_from_fault_params(L=params['L_km'], 
#                                                   z=params['z_km'], 
#                                                   dip=mc_d['dip_samp'], 
#                                                   D=max_eq_slip)
#
#    mc_d['max_M'] = eqs.calc_M_from_Mo(mc_d['max_Mo'])
#    
#    rec_int_df = pd.DataFrame(columns = rec_int_bins[1:],
#                              index=df_multi_ind, dtype='float')
#    t0 = time.time()
#    #print fdf.head(1),rec_int_df.head(1)
#    #print fdf.as_matrix().shape,len(mc_index)
#    print mc_index
#    out_list = Parallel(n_jobs=-3)( delayed( calc_iter_probs)(ii) 
#                                    for ii in mc_index)
#    print 'done with', fault, 'parallel calcs in {} s'.format((time.time()-t0))
#    for ii in mc_index:
#        fdf.loc[ii][:] = out_list[ii][0]
#        rec_int_df.loc[ii][:] = out_list[ii][1]
        
#    fdf.to_csv('../results/{}_test.csv'.format(fault))


# ### Calculate $P(M,t)$ for faults with characteristic frequency-magnitude distributions

# In[13]:


def calc_iter_probs(ii):
    df_iter = fdf.loc[ii].copy()
    df_iter['dip'] = mc_d['dip_samp'][ii]
    df_iter['Ddot'] = mc_d['Ddot_samp'][ii]

    # Generate EQ sample/sequence from F(M) dist.
    m_vec = np.linspace(5, mc_d['max_M'][ii], num=1000)
    fm_vec = eqs.F_char(m_vec, Mc=Mc, char_M=mc_d['char_M'][ii])
    M_samp = eqs.sample_from_pdf(m_vec, fm_vec, n_eq_samp)
    Mo_samp = eqs.calc_Mo_from_M(M_samp)
    
    # Make time series of earthquakes, including no eq years
    recur_int = eqs.calc_recurrence_interval(Mo=Mo_samp, 
                                             dip=mc_d['dip_samp'][ii],
                                             slip_rate=mc_d['Ddot_samp'][ii],
                                             L=params['L_km'],
                                             z=params['z_km'])

    cum_yrs = eqs.calc_cumulative_yrs(recur_int)
    eq_series = eqs.make_eq_time_series(M_samp, cum_yrs)
    
    # calculate probability of observing EQ in time_window
    for t in time_window:
        roll_max = pd.rolling_max(eq_series, t)
        df_iter[t] = (eqs.get_probability_above_value(roll_max, min_M_list)
                      * mc_d['dip_frac'])
        
    # calculate histgrams of recurrence intervals
    rec_int_counts_df = rec_int_df.loc[ii].copy()
    for mm in np.array(min_M_list):
        ints = np.diff( np.where(eq_series >= mm) )
        rec_int_counts_df.loc[mm] = np.histogram(ints, bins=rec_int_bins)[0]

    return df_iter,rec_int_counts_df

# In[ ]:

pl = False
if pl: 
    for fault in list(f.index):
        fdf = pd.DataFrame(index=df_multi_ind, columns=mc_cols, dtype='float')
        params = f.loc[fault]
        mc_d = {}
        mc_d['dip_samp'], mc_d['dip_frac'] = eqs.dip_rand_samp( params['dip_deg'], 
                                                             params['dip_err_deg'], 
                                                             mc_iters)

        mc_d['Ddot_samp'] = eqs.Ddot_rand_samp(params['slip_rate_mm_a'],
                                               params['sr_err_mm_a'], mc_iters)

        mc_d['max_Mo'] = eqs.calc_Mo_from_fault_params(L=params['L_km'], 
                                                       z=params['z_km'], 
                                                       dip=mc_d['dip_samp'], 
                                                       D=max_eq_slip)

        mc_d['max_M'] = eqs.calc_M_from_Mo(mc_d['max_Mo'])
        
        mc_d['char_Mo'] = eqs.calc_Mo_from_fault_params(L=params['L_km'], 
                                                        z=params['z_km'], 
                                                        dip=mc_d['dip_samp'], 
                                                        D=char_eq_slip)
    
        mc_d['char_M'] = eqs.calc_M_from_Mo(mc_d['char_Mo'])
        
        rec_int_df = pd.DataFrame(columns = rec_int_bins[1:],
                                  index=df_multi_ind, dtype='float')
        t0 = time.time()
        out_list = Parallel(n_jobs=n_cores)( delayed( calc_iter_probs)(ii) 
                                        for ii in mc_index)
        print 'done with', fault, 'parallel calcs in {} s'.format((time.time()-t0))
        for ii in mc_index:
            fdf.loc[ii][:] = out_list[ii][0]
            rec_int_df.loc[ii][:] = out_list[ii][1]
        
        fdf.to_csv('../results/{}_char.csv'.format(fault))


# ### Examining individual fault results

# #### Load datasets into Pandas dataframes
# 
# test with one:

# In[14]:


pv = pd.read_csv('../results/panamint_valley_gr.csv', index_col=[0])
pv_c = pd.read_csv('../results/panamint_valley_char.csv', index_col=[0])


# In[11]:


tw_xarray = np.tile(time_window, (mc_iters,1))

tw_xarray.shape


# In[12]:


tw_cols = list(time_window.astype('str'))
print tw_cols


# make some plotting functions

# In[21]:


def plot_probs(df, lw=0.5, alpha=0.2):

    df_5 = df[df.M==5]
    df_55 = df[df.M==5.5]
    df_6 = df[df.M==6]
    df_65 = df[df.M==6.5]
    df_7 = df[df.M==7]
    df_75 = df[df.M==7.5]

    aa = plot(tw_xarray.T, df_5[tw_cols].T, 'red', 
              lw=lw, alpha=alpha,label=r"$M\geq 5$")
    bb = plot(tw_xarray.T, df_55[tw_cols].T, 'green', 
              lw=lw, alpha=alpha,label=r"$M\geq 5.5$")
    cc = plot(tw_xarray.T, df_6[tw_cols].T, 'blue', 
              lw=lw, alpha=alpha,label=r"$M\geq 6$")
    dd = plot(tw_xarray.T, df_65[tw_cols].T, 'black', 
              lw=lw, alpha=alpha,label=r"$M\geq 6.5$")
    ee = plot(tw_xarray.T, df_7[tw_cols].T, 'gold', 
              lw=lw, alpha=alpha,label=r"$M\geq 7$")
    ff = plot(tw_xarray.T, df_75[tw_cols].T, 'orange', 
              lw=lw, alpha=alpha,label=r"$M\geq 7.5$")
    xlim(0,100)
    xlabel(r"$Time Window(t)$")
    ylabel(r"$P(M,t)$")

    return aa, bb, cc, dd, ee, ff


# In[19]:


def pmf_x_secs(df, tw, n_bins=101, hist_type='stepfilled', alpha=0.5):
    hist_bins = np.linspace(0, 1, num=n_bins)
    tw = str(tw)
    
    aa = df[df.M==5][tw].hist(bins=hist_bins, histtype=hist_type, 
                              color='r', alpha=alpha, label=r"$M\geq 5$")
    bb = df[df.M==5.5][tw].hist(bins=hist_bins, histtype=hist_type, 
                                color='g', alpha=alpha, label=r"$M\geq 5.5$")
    cc = df[df.M==6][tw].hist(bins=hist_bins, histtype=hist_type, 
                              color='b', alpha=alpha, label=r"$M\geq 6$")
    dd = df[df.M==6.5][tw].hist(bins=hist_bins, histtype=hist_type, 
                                color='k', alpha=alpha, label=r"$M\geq 6.5$")
    ee = df[df.M==7][tw].hist(bins=hist_bins, histtype=hist_type, 
                              color='gold', alpha=alpha, label=r"$M\geq 7$")
    ff = df[df.M==7.5][tw].hist(bins=hist_bins, histtype=hist_type, 
                                color='orange', alpha=alpha, label=r"$M\geq 7.5$")
    ylabel(r"Count")
    xlabel(r"$P(M,t=35)$")

    return aa, bb, cc, dd, ee, ff


# In[16]:


#plot_probs(pv)

#xlabel('observation time, yrs', fontsize=18)
#ylabel('probability', fontsize=18)
#title('Figure 1: Likelihood of observing an earthquake \n'
#      + 'larger than M on the Kongur Shan fault', fontsize=21)
#show()


# In[35]:


tw=35
nbins=1001
x_sec_lims = [0,0.4,0,400]

f = figure(dpi=120,figsize=(10,9))

subplot2grid((3,2), (0,0), colspan=2)
axvline(tw, lw=0.5, color='grey')
plot_probs(pv, lw=0.2, alpha=0.2)

subplot2grid((3,2), (1,0), colspan=2)
axvline(tw, lw=0.5, color='grey')
plot_probs(pv_c, lw=0.2, alpha=0.2)

subplot2grid((3,2), (2,0))
pmf_x_secs(pv, tw, n_bins=nbins)
axis(x_sec_lims)

subplot2grid((3,2), (2,1))
pmf_x_secs(pv_c, tw, n_bins=nbins)
axis([0,0.1,0,400])
f.subplots_adjust(wspace=0.5,hspace=0.5)
savefig("wprc_1.png")
close()
#show()


# The results for faults with a GR frequency-magnitude distribution
# indicate that it is unlikely that any individual fault would have an earthquake
# greater than $M \, 5$ in any modeled observation time (up to 100 years).  
# As an example, the results for the Panamint Valley fault are shown in
# Figure~\ref{fig:pv}a;
# this fault has the highest $P(M,t)$ of any of the well-studied LANFs.  The
# probability of observing a $\ge M \, 6.0$ event on the Panamint Valley fault 
# is about 0.5 for $t$ = 100 years, and about 0.15 for $t$ = 35 years, which is
# the length of the Global CMT catalog. As expected given the GR distribution,
# the $P(M,t)$ is much higher for smaller, more frequent events than for larger
# events.  
# 
# The results for faults with a characteristic frequency-magnitude distribution
# yield much lower $P(M,t)$ for small to moderate events, but for large events,
# $P(M,t)$ is higher (Figure~\ref{fig:pv}b,d); this is because the earthquake
# sequences are dominated by large, infrequent events, so the inter-event times
# for moderate events are several times greater. For the Panamint Valley fault,
# $P(M,t)$ for $ \ge M \, 5$ is about 0.07 (versus 0.25 for the GR distribution),
# but $P(M\ge 7, t=35)$ is around 0.025 (versus essentially zero for the GR
# distribution).  As the characteristic distribution likely better represents
# earthquakes on an individual large fault, these results suggest that is 
# very unlikely that we would expect to capture any significant seismicity
# on an single LANF in the moment tensor catalogs. A similar conclusion was 
# found by \citet{wernicke1995seis} based on a simple calculation, 
# assuming perfectly repeating large earthquakes on an idealized fault. 

# ### Earthquake Likelihood on All LANFs
# 
# To calculate the probability of observing an earthquake over the time
# window on any of the LANFs studied, we first assume that seismicity on
# each fault is independent and uncorrelated with seismicity on all
# other faults. This assumption is likely true for most faults, but may
# not be true for some proximal faults (such as the North and South
# Lunggar detachments, or the Papangeo and Tokorondo detachments),
# though it is unclear how to determine how these faults may interact 
# such that an appropriate conditional or joint probability may be calculated. 
# We determine the probability for each time window and minimum magnitude
# with the equation
# 
# \begin{equation}
# P_{AT \, or \, LP\, \ldots \, or \, DV} = 1 - (Q_{AT} \cdot Q_{LP} \cdot \ldots \, \cdot Q_{DV})
# %\label{ProbUnion}
# \end{equation}
# 
# where $P_{AT}$ is the probability of observing an earthquake on a
# single LANF (e.g., the Alto-Tiberina fault), and $Q_{AT} = 1 -
# P_{AT}$. Equation (\ref{ProbUnion}) is the union of probabilities for
# non mutually exclusive random events.

# ### Calculating these probabilities

# ####now load the rest:
# 
# Gutenberg-Richter models

# In[36]:


ks = pd.read_csv('../results/kongur_shan_gr.csv', index_col=0)
gm = pd.read_csv('../results/gurla_mandhata_gr.csv', index_col=0)
slr = pd.read_csv('../results/s_lunggar_gr.csv', index_col=0)
nlr = pd.read_csv('../results/n_lunggar_gr.csv', index_col=0)
pqxn = pd.read_csv('../results/pqx_n_gr.csv', index_col=0)
pqxq = pd.read_csv('../results/pqx_qingdu_gr.csv', index_col=0)
lp = pd.read_csv('../results/leo_pargil_gr.csv', index_col=0)
nqtl = pd.read_csv('../results/nqtl_gr.csv', index_col=0)
at = pd.read_csv('../results/alto_tiberina_gr.csv', index_col=0)
dv = pd.read_csv('../results/death_valley_gr.csv', index_col=0)
cd = pd.read_csv('../results/canada_david_gr.csv', index_col=0)
sd = pd.read_csv('../results/sevier_desert_gr.csv', index_col=0)
dxv = pd.read_csv('../results/dixie_valley_gr.csv', index_col=0)
dd = pd.read_csv('../results/dayman_dome_gr.csv', index_col=0)
pp = pd.read_csv('../results/papangeo_gr.csv', index_col=0)
tk = pd.read_csv('../results/tokorondo_gr.csv', index_col=0)
cb = pd.read_csv('../results/cordillera_blanca_gr.csv', index_col=0)
kz = pd.read_csv('../results/kuzey_gr.csv', index_col=0)
gn = pd.read_csv('../results/guney_gr.csv', index_col=0)


# characteristic models

# In[41]:


ks_c = pd.read_csv('../results/kongur_shan_char.csv', index_col=0)
gm_c = pd.read_csv('../results/gurla_mandhata_char.csv', index_col=0)
slr_c = pd.read_csv('../results/s_lunggar_char.csv', index_col=0)
nlr_c = pd.read_csv('../results/n_lunggar_char.csv', index_col=0)
pqxn_c = pd.read_csv('../results/pqx_n_char.csv', index_col=0)
pqxq_c = pd.read_csv('../results/pqx_qingdu_char.csv', index_col=0)
lp_c = pd.read_csv('../results/leo_pargil_char.csv', index_col=0)
nqtl_c = pd.read_csv('../results/nqtl_char.csv', index_col=0)
at_c = pd.read_csv('../results/alto_tiberina_char.csv', index_col=0)
dv_c = pd.read_csv('../results/death_valley_char.csv', index_col=0)
cd_c = pd.read_csv('../results/canada_david_char.csv', index_col=0)
sd_c = pd.read_csv('../results/sevier_desert_char.csv', index_col=0)
dxv_c = pd.read_csv('../results/dixie_valley_char.csv', index_col=0)
dd_c = pd.read_csv('../results/dayman_dome_char.csv', index_col=0)
pp_c = pd.read_csv('../results/papangeo_char.csv', index_col=0)
tk_c = pd.read_csv('../results/tokorondo_char.csv', index_col=0)
cb_c = pd.read_csv('../results/cordillera_blanca_char.csv', index_col=0)
kz_c = pd.read_csv('../results/kuzey_char.csv', index_col=0)
gn_c = pd.read_csv('../results/guney_char.csv', index_col=0)


# make list of faults

# In[37]:


f_list = ['ks', 'lp', 'gm', 'slr', 'nlr', 'pqxn', 'pqxq', 'nqtl', 'at',
          'dv', 'pv', 'cd', 'sd', 'dxv', 'dd', 'pp', 'tk', 'cb', 'kz', 'gn']

fault_list = [ks, lp, gm, slr, nlr, pqxn, pqxq, nqtl, at,
              dv, pv, cd, sd,  dxv, dd, pp, tk, cb, kz, gn]


# In[42]:


fault_c_list = [ks_c, lp_c, gm_c, slr_c, nlr_c, pqxn_c, pqxq_c, nqtl_c, at_c,
              dv_c, pv_c, cd_c, sd_c, dxv_c, dd_c, pp_c, tk_c, cb_c, kz_c, gn_c]

fc_list = ['ks_c', 'lp_c', 'gm_c', 'slr_c', 'nlr_c', 'pqxn_c', 
                'pqxq_c', 'nqtl_c', 'at_c', 'dv_c', 'pv_c', 'cd_c', 
                'sd_c', 'dxv_c', 'dd_c', 'pp_c', 'tk_c', 'cb_c', 'kz_c', 'gn_c']


# #### Calculate $q$, probability of *not* observing an earthquake, for each fault

# In[38]:


q = {}

for i, f_ in enumerate(f_list):
    q[f_] = fault_list[i].copy()
    q[f_][tw_cols] = 1 - q[f_][tw_cols]


# In[43]:


qc = {}

for i, f_ in enumerate(fc_list):
    qc[f_] = fault_c_list[i].copy()
    qc[f_][tw_cols] = 1 - qc[f_][tw_cols]


# #### Now estimate joint probabilities
# 
# Make list of columns to retain in final dataframe

# In[39]:


prob_cols = list( chain( *['M', list(tw_cols)] ) )
print prob_cols

# calculate $p_{joint}$ as $1-(q \cdot q \cdot q...)$ and make final dataframe

# In[40]:

#print q["pp"][tw_cols].as_matrix()
#all_prob_gr = 1 - np.product([qq[tw_cols] for qq in  q.values()])
#xx = np.array([qq[tw_cols].as_matrix() for qq in  q.values()])
#print xx.shape
#print type(qq[tw_cols])
all_prob_gr = pd.DataFrame()
#print q[f_list[0]]["5"]
L = len(q[f_list[0]]["5"])
print L
for tw in tw_cols:
    xx = np.ones(L)
    for f in f_list:
        xx = xx * np.array(q[f][tw])
    all_prob_gr[tw] = 1 - xx
#all_prob_gr = 1 - np.prod([qq[tw_cols].as_matrix() for qq in  q.values()])
#print all_prob_gr

all_prob_gr['M'] = ks['M'].tolist()

all_prob_gr = all_prob_gr.reindex_axis(prob_cols, axis=1)


# In[44]:
all_c_prob = pd.DataFrame()
#print q[f_list[0]]["5"]
L = len(qc[fc_list[0]]["5"])
print L
for tw in tw_cols:
    xx = np.ones(L)
    for f in fc_list:
        xx = xx * np.array(qc[f][tw])
    all_c_prob[tw] = 1 - xx
#all_prob_gr = 1 - np.prod([qq[tw_cols].as_matrix() for qq in  q.values()])
#print all_c_prob

all_c_prob['M'] = ks['M'].tolist()

all_c_prob = all_c_prob.reindex_axis(prob_cols, axis=1)



#all_c_prob = 1- np.product([qq[tw_cols] for qq in qc.values()])

#all_c_prob['M'] = ks['M']

#all_c_prob = all_c_prob.reindex_axis(prob_cols, axis=1)


# In[48]:


tw=35
nbins=501
x_sec_lims = [0,0.4,0,400]

f=figure(dpi=120,figsize=(10,9))

subplot2grid((3,2), (0,0), colspan=2)
axvline(tw, lw=0.5, color='grey')
plot_probs(all_prob_gr, lw=0.1, alpha=0.1)

subplot2grid((3,2), (1,0), colspan=2)
axvline(tw, lw=0.5, color='grey')
plot_probs(all_c_prob, lw=0.1, alpha=0.1)

subplot2grid((3,2), (2,0))
pmf_x_secs(all_prob_gr, tw, n_bins=nbins)
#axis(x_sec_lims)

subplot2grid((3,2), (2,1))
pmf_x_secs(all_c_prob, tw, n_bins=nbins)
#axis([0,0.1,0,400])
xlim([0,0.8])
f.subplots_adjust(wspace=0.5,hspace=0.5)
savefig("wprc_2.png")
close()
#show()


# In[50]:


tw=35
nbins=501
x_sec_lims = [0,0.4,0,400]

f=figure(dpi=120,figsize=(10,9))

subplot2grid((2,2), (0,0))
axvline(tw, lw=0.5, color='grey')
plot_probs(all_prob_gr, lw=0.1, alpha=0.1)

subplot2grid((2,2), (0,1))
axvline(tw, lw=0.5, color='grey')
plot_probs(all_c_prob, lw=0.1, alpha=0.1)

subplot2grid((2,2), (1,0))
pmf_x_secs(all_prob_gr, tw, n_bins=nbins)
#axis(x_sec_lims)

subplot2grid((2,2), (1,1))
pmf_x_secs(all_c_prob, tw, n_bins=nbins)
#axis([0,0.1,0,400])
xlim([0,0.8])
f.subplots_adjust(wspace=0.5,hspace=0.5)
savefig("wprc_3.png")
close()
#show()

