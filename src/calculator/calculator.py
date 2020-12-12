import random
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.stats import truncnorm
import json
from IPython.display import Image
import pandas as pd
from scipy.integrate import quad

#Used to generate values with normal distribution
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
  
def get_air_changes_per_hour(cfm, room_volume):
  if str(cfm) == 'nan':
      #impute unknown cfm with arbitrary cfm
      print('VAV unknown. Imputed with arbitrary VAV of 800 CFM')
      cfm = 800
  return (cfm * 60) / room_volume

#In the future generating occupants may be useful in simulating a indoor environment more accurately
def generate_occupants(n_occupants, n_infected, speaker_to_breathing_ratio, masks_efficacy = mask_efficacy):
    occupants = {}
    #Age input. If not age is inputted an age distribution should be generated for the region. 
    #Default age distribution is an empirical age distribution for a UCSD Classroom. 
    #Source: https://ir.ucsd.edu/_files/stats-data/profile/profile-2018-2019.pdf
    occupant_age = np.random.choice([18, 19, 20, 21, 22, 23, 24, random.randint(25, 30)], size = n_occupants, p = [.198, .195, .195, .1875, .09375, .046875, .023875, .06])
    for i in range(n_occupants):
        curr_occ = "occ" + str(i)
        occupants[curr_occ] = {}
        occupants[curr_occ]['age'] = occupant_age[i]
        #Here we must change respiratory function to reflect exercise
        #Normal breathing/speaking/singing concentrations
        #Concentration from breathing Standard: 0.1 Range: 0.06–1.0 cubic centimeters Distribution: Normal
        #Source:
        normal_breathing = get_truncated_normal(mean=.1, sd=.02, low=.06, upp=1)
        occupants[curr_occ]['conc_breathing'] = normal_breathing.rvs() * CUBIC_CM_TO_METERS
        #Concentration from speaking (singing) Standard: 1.1 Range: 0.06–6.0 cubic centimeters Distribution: Normal
        #Source:
        normal_speaking = get_truncated_normal(mean=1.1, sd=.02, low=.06, upp=6)
        occupants[curr_occ]['conc_speaking'] = normal_speaking.rvs() * CUBIC_CM_TO_METERS
        #Respiratory rate standard: 10 Range: 5–20 L/min Distribution: Normal
        #Source:
        normal_respiratory_rate = get_truncated_normal(mean= 10, sd=5, low=5, upp=20)
        occupants[curr_occ]['norm_resp_rate'] = normal_respiratory_rate.rvs()
        occupants[curr_occ]['mask_efficacy'] = masks_efficacy
        occupants[curr_occ]['testedWeekly'] = True
    return occupants
  
  
#Calculate emission rates
def get_quanta_emmission_rate(cv, ci, IR, Dc, Dv = droplet_vol):
    #Convert droplet volume from cubic micrometers to centimeters
    summation = sum([Dc['.8μm'] * (Dv['.8μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['1.8μm'] * (Dv['1.8μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['3.5μm'] * (Dv['3.5μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['5.5μm'] * (Dv['5.5μm'] * CUBIC_μM_TO_CUBIC_CM)])
    #Convert IR from cubic meters to mililiter
    return cv * ci * (IR * CUBIC_M_TO_ML) * summation
  
def quanta_concentration(t, I = n_infected, ERq = ERq, IVRR = ivrr, V = room_volume_m, n0 = initial_quanta):
  return ((ERq * I) / (IVRR * V)) + (n0 + ((ERq * I) / IVRR)) * ((np.e**(-IVRR * t)) / V)

def calculate_risk(t, IR = inhalation_rate[activity]):
  ans, err = quad(quanta_concentration, 0, t)
  return 1 - np.e**(-IR * ans)

#Calculate maximum people allowed in the room given an exposure time (hours)
#steady state model
def calc_n_max_ss(exp_time, relative_humidity, max_aerol_radius,max_viral_deact_rate, room_vol, air_exch_rate, mean_ceiling_height_m, ERq, mask): #exp time in hrs
    eff_aerosol_radius = ((0.4 / (1 - relative_humidity)) ** (1 / 3)) * max_aerosol_radius
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * relative_humidity
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 5
    air_filt_rate = aerosol_filtration_eff[1] * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    n_max = 1 + 0.1 / (airb_trans_rate * exp_time)
    return n_max

#transient model
def calc_n_max_t(exp_time, relative_humidity, max_aerol_radius,max_viral_deact_rate, room_vol, air_exch_rate, mean_ceiling_height_m, ERq,mask): #exp time in hrs
    eff_aerosol_radius = ((0.4 / (1 - relative_humidity)) ** (1 / 3)) * max_aerosol_radius
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * relative_humidity
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 5
    air_filt_rate = aerosol_filtration_eff[1] * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    n_max = 1 + (0.1 * (1 + 1/(conc_relax_rate * exp_time)) / (airb_trans_rate * exp_time))
    return n_max

#Calculate maximum exposure time allowed given a capacity (# people):
def calc_max_time(n_max, relative_humidity, max_aerol_radius,max_viral_deact_rate, room_vol, air_exch_rate, mean_ceiling_height_m, ERq,mask):
    #risk_tolerance = risk_tolerance
    eff_aerosol_radius = ((0.4 / (1 - relative_humidity)) ** (1 / 3)) * max_aerosol_radius
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * relative_humidity
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 5
    air_filt_rate = aerosol_filtration_eff[1] * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    exp_time_ss = 0.1 / ((n_max - 1) * airb_trans_rate)  # hrs, steady-state
    exp_time_trans = exp_time_ss * (1 + (1 + 4 / (conc_relax_rate * exp_time_ss)) ** 0.5) / 2  # hrs, transient
    return exp_time_trans
