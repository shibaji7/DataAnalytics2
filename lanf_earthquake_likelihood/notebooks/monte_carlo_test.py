
# coding: utf-8

# # Monte Carlo simulation to find recurrence interval PDF on a fault

# In[1]:


#get_ipython().magic(u'pylab inline')
#get_ipython().magic(u"config InlineBackend.figure_format = 'retina'")


# In[2]:

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# In[3]:


f = pd.read_csv('../data/lanf_stats.csv', index_col=0)


# ## Define and sample from frequency-magnitude distribution F(M)

# ### Define F(M)

# In[4]:


def F(M, Mmin=2.5, Mc=7.64, B=0.65):
    """
    F(M) is the tapered Gutenberg-Richter distribution
    given by Shen et al., 2007 SRL with values for constants
    from the CRB values of Bird and Kagan 2004 BSSA.
    """
    term1 = 10**(-1.5 * B * (M - Mmin) )
    term2 = np.exp(10**(1.5 * (Mmin - Mc) - 1.5 * (M-Mc) ) )
    
    return term1 * term2


# In[5]:


M_vec = np.linspace(5, 7.64, num=1000)


# In[6]:


F_Mvec = F(M_vec)


# In[7]:


plt.plot(M_vec, F_Mvec)
plt.xlabel('M (vals)')
plt.ylabel('Freq (counts)')
plt.savefig("mcmc_1.png")

# ### Sample from F(M) with inverse transform sampling algorithm

# In[8]:

plt.close()

plt.plot(M_vec, np.cumsum(F_Mvec) )
plt.xlabel('M (vals)')
plt.ylabel('Cumulative Freq (counts)')
plt.savefig("mcmc_2.png")
plt.close()

# #### Define an algorithm to do inverse transform sampling.
# 
# In order to not have to fit a curve to the ECDF for sampling,
# we will instead round incoming random values to the points
# closest to them in the ECDF.  Shady, I know, but...  It also
# allows for sampling from arbitrary PDFs, not ones that just look
# like some simple function.

# In[9]:


def get_ecdf(counts, norm=True):
    """Calculates an array sort of like an ECDF, but no x reference"""
    if norm == True:
        ecdf = np.cumsum(counts) / np.float(np.sum(counts))
    else:
        ecdf = np.cumsum(counts)
    
    return ecdf


def get_ecdf_dict(vals=None, counts=None, norm=True):
    """
    Returns a dict object linking counts to values.
    Things better be sorted!
    """
    ecdf = get_ecdf(counts, norm=norm)
    ecdf_dict = {}
    
    for i, val in enumerate(vals):
        ecdf_dict[val] = ecdf[i]
    
    return ecdf_dict


def get_inverse_ecdf_dict(vals=None, counts=None, norm=True):
    
    ecdf_dict = get_ecdf_dict(vals, counts, norm)
    inv_ecdf_dict = {y:x for x,y in ecdf_dict.iteritems()}
    
    return inv_ecdf_dict
    

def find_nearest_dict_val(dict_, vals):
    key_array = np.array(dict_.keys())
    nearest_vals = np.empty(vals.shape)
    
    for i, val in enumerate(vals):
        idx = (np.abs(key_array - val)).argmin()
        nearest_vals[i] = dict_[ key_array[idx] ]
    
    return nearest_vals


# In[10]:


plt.plot(M_vec, get_ecdf(F_Mvec))
plt.xlabel('M (vals)')
plt.ylabel('ECDF Freq (counts)')
plt.savefig("mcmc_2.png")
plt.close()

# In[11]:


ECDF_dict = get_ecdf_dict(vals=M_vec, counts=F_Mvec, norm=True)


# In[12]:


plt.plot(ECDF_dict.keys(), ECDF_dict.values(), ',')
plt.xlabel('M (vals)')
plt.ylabel('Cumulative Freq (counts)')
plt.savefig("mcmc_2.png")


# In[13]:


inv_ecdf_dict = get_inverse_ecdf_dict(vals=M_vec, counts=F_Mvec, norm=True)


# In[14]:


rand_samp = np.random.rand(int(1e4))


# In[15]:


freq_mag_samp = find_nearest_dict_val(inv_ecdf_dict, rand_samp)


# In[16]:

plt.close()
plt.hist(freq_mag_samp, bins=100)
plt.xlabel('Mag samples')
plt.ylabel('Count')
plt.xlim(5,7.5)
plt.savefig("mcmc_3.png")
plt.close()

# ## Trial for South Lunggar

# #### Get info for South Lunggar Detachment

# In[16]:


slr = f.loc['s_lunggar']


# In[17]:


slr


# #### Make dataframe for S. Lunggar MC trial

# In[18]:


n_rand = int(1e5)


# In[19]:


SLR_MC = pd.DataFrame(index=np.arange(n_rand), columns=['M', 'Mo', 'Ddot', 'dip', 'recur_int'])


# #### Get random EQ values (M and Mo) from F(M) distribution

# In[20]:


FM_icdf_d = get_inverse_ecdf_dict(vals=M_vec, counts=F_Mvec, norm=True)

SLR_MC['M']= find_nearest_dict_val(FM_icdf_d, np.random.rand(n_rand) )


# In[21]:


def Mo_from_M(M, C=6):
    """
    Calculate seismic moment (Mo) from
    moment magnitude (M) given a scaling law
    """
    term1 = 3/2. * C * (np.log(2) + np.log(5) )
    term2 = 3/2. * M * (np.log(2) + np.log(5) )
    
    Mo = np.exp( term1 + term2)
    
    return Mo


# In[22]:


SLR_MC['Mo'] = Mo_from_M(SLR_MC.M)


# #### Get random samples from slip rate and dip distributions

# In[23]:


SLR_MC['Ddot'] = (np.random.rand(n_rand) * 2 * slr['sr_err_mm_a'] 
                 + slr['slip_rate_mm_a'] - slr['sr_err_mm_a'] )


# In[25]:


plt.hist(SLR_MC.Ddot, bins=100)
plt.xlabel("Ddot")
plt.savefig("out/mcmc_4.png")
plt.close()

# In[24]:


SLR_MC.dip = (np.random.rand(n_rand) * 2 * slr['dip_err_deg'] 
             - slr['dip_err_deg'] + slr['dip_deg'] )


# #### Calculate recurrence interval from Mo, fault geometry, slip rate

# In[25]:


def calc_rec_int(Mo=None, dip=None, mu=6e9, L=None, z=None,
                 slip_rate=None):
    
    return Mo * np.sin(dip) / (mu * L * z * slip_rate)


# In[26]:


SLR_MC.recur_int = calc_rec_int(Mo=SLR_MC.Mo, dip=np.radians(SLR_MC.dip), 
                                L=slr['L_km']*1000, z=slr['z_km']*1000, 
                                slip_rate=SLR_MC.Ddot/1000.)


# In[27]:


np.median(SLR_MC.recur_int)


# In[28]:


np.mean(SLR_MC.recur_int)


# In[29]:


np.median(SLR_MC.M)


# In[32]:


plt.plot(SLR_MC.M, SLR_MC.recur_int, '.')
plt.xlabel('M')
plt.ylabel('recurrence int (yr)')
plt.savefig("mcmc_5.png")
plt.close()
# #### Calculate mean slip (D) for each event

# In[30]:


SLR_MC['D_m'] = SLR_MC.Ddot / 1000. * SLR_MC.recur_int


# In[31]:


SLR_MC.head()


# In[35]:


SLR_MC['cum_yrs'] = np.int_(np.cumsum( SLR_MC.recur_int.round() ) )


# In[36]:


SLR_MC.head()


# In[37]:


SLR_MC.cum_yrs.max()


# from here I think I need to brute force the fuck out of things.
# 
# - Make a time series of EQs.  Zeros will be years with no EQs.  Non-zero years will be earthquake magnitudes.  The time from one earthquake to the next is the amount of time required to build up sufficient moment for the magnitude of the later earthquake.
# - Then do a rolling maximum with a given time window (indicating period of observation) to calculate the number of time windows that have events larger than some magnitude threshold.
# 
# Maybe there is a better way to go about it, but I don't understand the statistical treatments.

# In[38]:


eq_time_series = np.zeros(SLR_MC.cum_yrs.max())


# In[39]:


for i in np.arange(n_rand-1):
    yr = SLR_MC.cum_yrs[i]
    mag = SLR_MC.M[i]
    eq_time_series[yr] = mag


# In[40]:


eq_max_40yr = pd.rolling_max(eq_time_series, 40)


# In[41]:


def prob_above_val(series, val):
    count_above = len( series[series >= val])
    
    return count_above / float(len(series) )


# In[42]:


prob_above_val(eq_max_40yr, 5.5)


# In[43]:


window_prob_d = {}


# In[44]:


min_eq_size = 6.5


# In[45]:


for i in np.arange(100) + 1:
    eq_max_series = pd.rolling_max(eq_time_series, i)
    window_prob_d[i] = prob_above_val(eq_max_series, min_eq_size)


# In[46]:


plt.plot(window_prob_d.keys(), window_prob_d.values())
plt.xlabel('observation window (yrs)')
plt.ylabel('prob. observing eq >= {}'.format(min_eq_size) )
plt.savefig("mcmc_6.png")
plt.close()

# I think this is works.

# ## Define new F(M) that includes characteristic events

# In[32]:


def lognormal(x, mu = -0.5, sigma=0.5):
    term1 = 1 / sigma * x * np.sqrt(2 * np.pi)
    term2 = np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2) )
    
    return term1 * term2


def F_char(M, Mmin=2.5, Mc=7.64, B=0.65, char_M=6.25,
           char_amplitude_scale=5.25):
    """
    F(M) is the tapered Gutenberg-Richter distribution
    given by Shen et al., 2007 SRL with values for constants
    from the CRB values of Bird and Kagan 2004 BSSA.
    """
    Fm = F(M, Mmin=Mmin, Mc=Mc, B=B)
    F_char_amp = F(char_amplitude_scale)
    
    M_shift = M - char_M + 1
    M_shift[M_shift < 0] = 0
    
    F_lognormal = lognormal( (M_shift) )
    norm_constant = np.max(F_lognormal) / F_char_amp
    
    F_lognormal *= 1/ norm_constant
    
    F_char = Fm + F_lognormal
    
    F_char += - F_char.min()
    
    return F_char



# In[33]:


plt.plot(M_vec, F(M_vec) ,label="GR")


# In[34]:


plt.plot(M_vec, F_char(M_vec), label="Char")
plt.xlabel("Mag")
plt.ylabel("Relative Frequency")
plt.xlim(5,7.5)
plt.savefig("mcmc_7.png")
plt.close()


# schweet.

# ### Now let's do this for S. Lunggar again.

# Make GB+characteristic F(M) distribution, and sample from it

# In[36]:


SLR_CH_MC = pd.DataFrame(index=np.arange(n_rand), columns=SLR_MC.columns)


# In[37]:


FM_char = F_char(M_vec, char_amplitude_scale=5.25)


# In[38]:


FM_char_icdf_d = get_inverse_ecdf_dict(vals=M_vec, counts=FM_char, norm=True)

SLR_CH_MC.M = find_nearest_dict_val(FM_char_icdf_d, np.random.rand(n_rand) )


# In[39]:


plt.plot(M_vec, FM_char)
plt.legend(loc=0)
plt.xlabel("Mag")
plt.ylabel("Relative Frequency")
plt.xlim(5,7.5)
plt.savefig("mcmc_71.png")
plt.close()

# In[40]:


plt.plot(M_vec, get_ecdf(FM_char) , 'b')
plt.savefig("mcmc_8.png")
plt.close()

plt.plot(M_vec, get_ecdf(F_Mvec), 'r')
plt.savefig("mcmc_9.png")
plt.close()

# In[41]:


SLR_CH_MC.M.hist(bins=100)


# Calculate Mo and recurrence interval distribution

# In[73]:


SLR_CH_MC.Mo = Mo_from_M(SLR_CH_MC.M)

SLR_CH_MC['Ddot'] = (np.random.rand(n_rand) * slr['sr_err_mm_a'] 
                     - slr['sr_err_mm_a']/2. + slr['slip_rate_mm_a'] )

SLR_CH_MC.dip = (np.random.rand(n_rand) * slr['dip_err_deg'] 
                 - slr['dip_err_deg']/2. + slr['dip_deg'] )


# In[74]:


SLR_CH_MC.recur_int = calc_rec_int(Mo=SLR_CH_MC.Mo, dip=np.radians(SLR_CH_MC.dip), 
                                  L=slr['L_km']*1000, z=slr['z_km']*1000, 
                                  slip_rate=SLR_CH_MC.Ddot/1000.)


# In[75]:


np.median(SLR_CH_MC.recur_int)


# In[76]:


np.mean(SLR_CH_MC.recur_int)


# In[77]:


SLR_CH_MC['cum_yrs'] = np.int_(np.cumsum( SLR_CH_MC.recur_int.round() ) )


# In[78]:


SLR_CH_MC.cum_yrs.max()


# In[79]:


ch_eq_time_series = np.zeros(SLR_CH_MC.cum_yrs.max())


# In[80]:


for i in np.arange(n_rand-1):
    yr = SLR_CH_MC.cum_yrs[i]
    mag = SLR_CH_MC.M[i]
    ch_eq_time_series[yr] = mag


# In[81]:


ch_window_prob_d = {}


# In[82]:


for i in np.arange(100) + 1:
    eq_max_series = pd.rolling_max(ch_eq_time_series, i)
    ch_window_prob_d[i] = prob_above_val(eq_max_series, min_eq_size)


# In[83]:


plt.plot(ch_window_prob_d.keys(), ch_window_prob_d.values())
plt.xlabel('observation window (yrs)')
plt.ylabel('prob. observing eq >= {}'.format(min_eq_size) )
plt.title('Probability of observing EQ in time window, \n Characteristic F(M) dist')
plt.savefig("mcmc_10.png")
plt.close()

# In[84]:


plt.plot(ch_window_prob_d.keys(), ch_window_prob_d.values(), 'b', label='char')
plt.plot(window_prob_d.keys(), window_prob_d.values(), 'r', label='exp')
plt.xlabel('observation window (yrs)')
plt.ylabel('prob. observing eq >= {}'.format(min_eq_size) )
plt.legend(loc='upper left')
plt.savefig("mcmc_11.png")
plt.close()

#show()


# This is a result.  Is it intuitive, that we are less likely to see a large earthquake over some time series if the earthquakes have a characteristic F(M) distribution vs. an exponential (GB) distribution?
# 
# I was a little surprised, but I think it has to do with the way the sampling is done (which I would not say is wrong).  Basically, my interpretation is that the time series is longer in the characteristic sequence (by 2 million years), so there is more 'empty' time, meaning that the probability of seeing an earthquake is less.
# 
# But I don't know...
# 
# I think that the differences are small enough that it doesn't matter.

# In[60]:




