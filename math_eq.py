import math
import numpy as np
#formual to get distance in parsecs
def get_distance_in_parsecs(parallax):
    distance_psc = 1000 / parallax
    return distance_psc
#formula to get the absolute magnitude in visible wavelength    
def get_absolute_magnitude_visible(v_band, 
                                   get_distance_in_parsecs, 
                                   v_extinction):
    M_v = v_band - ((5 * math.log10(get_distance_in_parsecs) + 
                             v_extinction - 5))
    return M_v
#y = mx + c for period luminsoity relationship 
def get_period_luminosity_relation(P, a, b):
    return a*P + b
#formual to get abosolute magnitude in infrared
def get_absolute_magnitude_infra_red(i_band, 
                                     get_distance_in_parsecs,
                                     v_extinction):
    M_i = i_band - ((5 * math.log10(get_distance_in_parsecs) + 
                             v_extinction*0.556 - 5))
    return M_i
#formual to get log period from data    
def get_log_period(period):
    log_period = math.log10(period)
    return log_period
#from lab books (propogation of error formuula)
#z = 1/A, error = (error/A**2)
def get_distance_err(p_err, parallax):
    distance_err = 1000*(p_err/(parallax**2))
    return distance_err
#formula to get absolute magnitude in visible
def get_absolute_magnitude_visible_error(get_distance_err, 
                                         get_distance_in_parsecs, 
                                         av_err):
    err = math.sqrt(((5 * get_distance_err / (get_distance_in_parsecs * 
                                              math.log(10)))) ** 2 + (av_err)**2)
    return err
#formual to get absolute magnitude infrared
def get_absolute_magnitude_infra_red_error(get_distance_err, 
                                           get_distance_in_parsecs, 
                                           av_err):
    err = math.sqrt((5 * get_distance_err / (get_distance_in_parsecs * 
                                                       math.log(10)))**2 + (av_err*0.556)**2)
    
    return err
#getting distance of stars from a*P + b and apparent magnitude and extinction
def get_distance_visible(m_av, visible_alpha, log_p, visible_beta, a_v):
    distance_v = (10**((m_av-(visible_alpha * log_p + visible_beta) + 5 - a_v)/5))/10**6
    return distance_v
#same as above but in infrared
def get_distance_ir(m_air, ir_alpha, log_p, ir_beta, a_v):
    distance_ir = (10**((m_air-(ir_alpha * log_p + ir_beta) + 5 - (a_v*0.556))/5))/10**6
    return distance_ir
#ipythonnotebook lesson 13: formula for chi squared
def chisq(observed, expected, expected_i):
    chi = np.sum(((observed-expected)**2.0)/(expected_i**2.0))
    return chi
#formula to propogate the error in distance 
def error_propogation(x, sigma_m, sigma_b):
    sigma_sqrd = np.sqrt((x**2)*(sigma_m**2) + sigma_b**2)
    return sigma_sqrd
#spliting data for each each star by galaxy and adding it to a array total distance
def split_data(star_data, length):
    b = 0
    total_distance = []
    for i in range(len(length)):
        distance = []
        for a in range(length[i]):
            distance.append(star_data[b])
            b+=1
        total_distance.append(distance)
    return total_distance    
#formula of inverse variance weighted mean to aggregate the data 
#to minimize the variance of the weighted average
def inverse_mean(data, error):
    return np.sum(np.array(data)/np.array(error)**2)/np.sum(1/np.array(error)**2)
#formula to get the error
def inverse_mean_error(error):
    return np.sqrt(1/np.sum(1/np.array(error)**2))
#y = mx for hubble's law
def hubbleLaw(H0,d):
    return H0*d