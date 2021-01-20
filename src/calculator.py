import random
import numpy as np
import sys, os



import json

import pandas as pd
from scipy.integrate import quad

MINUTES_IN_HOUR = 60
SECONDS_IN_MINUTE = 60
SECONDS_IN_HOUR = 3600
FEET_TO_METERS = 0.3048
SQ_FT_TO_SQ_M = 0.092903
CUBIC_FT_TO_METERS = 0.0283168
CUBIC_CM_TO_METERS = 1e-6
CUBIC_μM_TO_CUBIC_CM = 1e-12
CUBIC_M_TO_ML = 1e6

var = {
"wet_asl_d": 5,
"virus_life": 1.7,
"hi_viral_load": 500000000,
"si_viral_load": 5000000000,
"deposition_prob": 0.5,
"cv": 1e9,
"ci": 0.02,
"IR": {"resting": 0.49,
       "standing": 0.54,
       "light_exercise": 1.38,
       "moderate_exercise": 2.35,
       "heavy_exercise": 3.30},
"droplet_conc": {
                "speaking": {".8μm": 0.4935,"1.8μm": 0.1035, "3.5μm": 0.073, "5.5μm": 0.035},
                "counting":  
                            {".8μm": 0.236, "1.8μm": 0.068, "3.5μm": 0.007, "5.5μm": 0.011},
                "whispering":
                            {".8μm": 0.110, "1.8μm": 0.014, "3.5μm": 0.004, "5.5μm": 0.002},
                "singing":
                            {".8μm": 0.751, "1.8μm": 0.139, "3.5μm": 0.139, "5.5μm": 0.059},
                "breathing":
                            {".8μm": 0.084, "1.8μm": 0.009, "3.5μm": 0.003, "5.5μm": 0.002}},
    
"droplet_vol": {".8μm": 0.26808257310632905, "1.8μm": 3.053628059289279, "3.5μm": 22.44929750377706, "5.5μm": 87.11374629016697},
"pass_vent_rate" : 0.35,
"deposition_rate" : 0.24,
    
    
"viral_inactivation" : 0.63,
"initial_quanta" : 0,
"viral_load_sptm" : 1000000000
}

#Used to generate values with normal distribution

#Helper Functions
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
    
def enablePrint():
    sys.stdout = sys.__stdout__

def get_air_changes_per_hour(cfm, room_volume):
    if str(cfm) == 'nan':
        #impute unknown cfm with arbitrary cfm
        print('VAV unknown. Imputed with arbitrary VAV of 800 CFM')
        cfm = 800
    elif cfm == 0:
        #Natural ACH of .15 if CFM is off provided
        return .15
    return (cfm * 60) / room_volume

def get_room_data(filepath, room_id):
    CUBIC_FT_TO_METERS = 0.0283168
    room_table = pd.read_csv(filepath)
    room_dic = {}
    if len(room_table.loc[room_table['Room'] == room_id]) == 0:
        return print("User input error: Room " + room_id + " not found.")
    room_table.loc[room_table['Room'] == room_id]['Area']
    #Room Area in square ft.
    room_dic['room_area'] = room_table.loc[room_table['Room'] == room_id]['Area'].item()
    #Room Height. Average room height of 10 ft is chosen if nan
    room_hght = room_table.loc[room_table['Room'] == room_id]['Height'].item()
    if np.isnan(room_hght):
        print(room_id + ' Room height not found. Average room height of 10 ft imputed.')
        room_hght = 10
    room_dic['room_hght'] = room_hght
    
    #CFM range. If no CFM is provided min is chosen by default
    cfm_range = room_table.loc[room_table['Room'] == room_id]['VAV'].item()
    if isinstance(cfm_range, float):
        ca_requirement_cfm = .53 * room_dic['room_area']
        room_dic['cfm_range'] = [ca_requirement_cfm, ca_requirement_cfm]
        print(room_id + ' CFM rate not found. California minimum ventilation requirement imputed')
    else:
        room_dic['cfm_range'] = list(map(int, cfm_range.split(',')))
    
    #Windows
    room_dic['windows'] = room_table.loc[room_table['Room'] == room_id]['Windows'].item()
    
    #V is volume of room
    room_dic['room_volume'] = room_dic['room_area'] * room_hght
    #Unit Conversion
    room_dic['room_volume_m'] = room_dic['room_area'] * room_hght * CUBIC_FT_TO_METERS
    
    return room_dic

def get_quanta_emmission_rate(activity, expiratory_activity, var = var):
    CUBIC_μM_TO_CUBIC_CM = 1e-12
    CUBIC_M_TO_ML = 1e6
    Dc = var['droplet_conc'][expiratory_activity]
    Dv = var['droplet_vol']
    #Convert droplet volume from cubic micrometers to centimeters
    summation = sum([Dc['.8μm'] * (Dv['.8μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['1.8μm'] * (Dv['1.8μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['3.5μm'] * (Dv['3.5μm'] * CUBIC_μM_TO_CUBIC_CM),
                     Dc['5.5μm'] * (Dv['5.5μm'] * CUBIC_μM_TO_CUBIC_CM)])
    #Convert IR from cubic meters to mililiter
    return var['cv'] * var['ci'] * (var['IR'][activity] * CUBIC_M_TO_ML) * summation
      
#Infection Risk Calculator
def infection_risk(t, room_id, n_occupants, activity, expiratory_activity, room_data_path, var = var):
    #CFM can be boolean or number within range of cfm_range
    CUBIC_μM_TO_CUBIC_CM = 1e-12
    ERq = get_quanta_emmission_rate(activity, expiratory_activity)
    room_dic = get_room_data(room_data_path, room_id)
    cfm_range = room_dic['cfm_range']
    cfm = sum(cfm_range) / 2
    #Air Changes per Hour
    air_change_rate = get_air_changes_per_hour(cfm, room_dic['room_volume'])
    
    ##To calculate infection rate we will aggregate the past week of testing for UC San Diego (last updated: 12/10/20)
    #Source: https://returntolearn.ucsd.edu/dashboard/index.html
    infection_rate = (2 + 11 + 10+ 3 + 7 + 10 + 12 + 5)/(20 + 1385 + 1375 + 286 + 1332 + 1414 + 944 + 1244)
    n_infected = infection_rate * n_occupants
    if n_infected < 1:
        n_infected = 1
    #Infectious virus removal rate
    ivrr = air_change_rate + var['deposition_rate'] + var['viral_inactivation']
    
    def quanta_concentration(t, I = n_infected, ERq = ERq, IVRR = ivrr, V = room_dic['room_volume_m'], n0 = var['initial_quanta']):
        return ((ERq * I) / (IVRR * V)) + (n0 + ((ERq * I) / IVRR)) * ((np.e**(-IVRR * t)) / V)

    ans, err = quad(quanta_concentration, 0, t)
    
    risk = 1 - np.e**(-var['IR'][activity] * ans)
    
    print('The resulting risk of infection is ' + str(risk * 100) +'%')
    print('It is predicted that ' + str(risk) + ' x ' + str(n_occupants) + ' = ' + str(int(risk * n_occupants)) + ' susceptible occupants will be infected')
    
    
    return risk




#Calculate maximum people allowed in the room given an exposure time (hours)
#steady state model
def calc_n_max_ss(exp_time, ax_aerosol_radius, room_area,room_height, air_exch_rate,IR,Dc,DV,mask): #exp time in hrs
    cv = 100
    ci = 1* 10 **(max(range(5, 10)))
    ERq = get_quanta_emmission_rate(cv, ci, IR, Dc, Dv)
    room_vol = room_height * room_area
    room_vol_m = room_vol*0.0283168
    mean_ceiling_height_m = room_height*0.3048
    eff_aerosol_radius = ((0.4 / (1 - 0.4)) ** (1 / 3)) * max_aerosol_radius
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * 0.4
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 10
    air_filt_rate = 0.1 * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    n_max = 1 + 0.1 / (airb_trans_rate * exp_time)
    return n_max

#transient model
def calc_n_max_t(exp_time,max_aerosol_radius, room_area,room_height, air_exch_rate,IR,Dc,DV,mask): #exp time in hrs
    cv = 100
    ci = 1* 10 **(max(range(5, 10)))
    ERq = get_quanta_emmission_rate(cv, ci, IR, Dc, Dv)
    eff_aerosol_radius = ((0.4 / (1 - 0.4)) ** (1 / 3)) * max_aerosol_radius
    room_vol = room_height * room_area
    room_vol_m = room_vol*0.0283168
    mean_ceiling_height_m = room_height*0.3048
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * 0.4
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 10
    air_filt_rate = 0.1 * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    n_max = 1 + (0.1 * (1 + 1/(conc_relax_rate * exp_time)) / (airb_trans_rate * exp_time))
    return n_max

#Calculate maximum exposure time allowed given a capacity (# people):
def calc_max_time(n_max, max_aerosol_radius, room_area,room_height, air_exch_rate,IR,Dc,DV,mask):
    cv = 100
    ci = 1* 10 **(max(range(5, 10)))
    ERq = get_quanta_emmission_rate(cv, ci, IR, Dc, Dv)
    room_vol = room_height * room_area
    room_vol_m = room_vol*0.0283168
    mean_ceiling_height_m = room_height*0.3048
    eff_aerosol_radius = ((0.4 / (1 - 0.4)) ** (1 / 3)) * max_aerosol_radius
    sett_speed_mm = 3 * (eff_aerosol_radius / 5) ** 2 #mm/s
    sett_speed = sett_speed_mm * 60 * 60 / 1000  # m/hr
    viral_deact_rate = 0.3 * 0.4
    fresh_rate = room_vol * air_exch_rate / 60
    recirc_rate = fresh_rate * (1/0.5 - 1)
    exhaled_air_inf = ERq * 10
    air_filt_rate = 0.1 * recirc_rate * 60 / room_vol #have to specify which filtration we have
    conc_relax_rate = air_exch_rate + air_filt_rate + viral_deact_rate + sett_speed / mean_ceiling_height_m
    airb_trans_rate = ((0.5 * mask) ** 2) * exhaled_air_inf / (room_vol_m * conc_relax_rate)
    exp_time_ss = 0.1 / ((n_max - 1) * airb_trans_rate)  # hrs, steady-state
    exp_time_trans = exp_time_ss * (1 + (1 + 4 / (conc_relax_rate * exp_time_ss)) ** 0.5) / 2  # hrs, transient
    return exp_time_trans
