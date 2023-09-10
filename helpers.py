#Calder Blackman
#created 7/28/23

#TO-DO tags denote incomplete routines
#UPDATE tags denote code that probably needs to be updated

#HCD tags denote hardcoded values/parameters

#math
import math
import numpy as np

#returns dist between two bodies (see class above); must pass type Body (or Observer actually)
def dist(b1, b2):
    return np.linalg.norm(np.array([b1.get_xyz()]) - np.array([b2.get_xyz()]))

#calc app mag given observer, band of interest, sun, and body
#expects observer to be type Observer and band to be a string
#expects sun and body to be type Body
def get_app_mag(observer, band, sun, body):
    #sun,body distance
    d_sb = dist(sun, body)
    #print(f"d_sb = {d_sb}")
    
    #body,observer distance
    d_bo = dist(body, observer)
    #print(f"d_bo = {d_bo}")

    #observer,sun distance
    d_os = dist(observer, sun)
    #print(f"d_os = {d_os}")
    
    #law of cosines
    phase_angle = math.acos((d_bo**2 + d_sb**2 - d_os**2) / (2 * d_bo * d_sb))
    #print(f"phase_angle = {phase_angle}")

    #phase function: assumes ideal diffuse reflecting sphere
    phase_factor = (2 / 3) * ((1 - phase_angle / math.pi) * math.cos(phase_angle) + 
                                (1 / math.pi) * math.sin(phase_angle))
    #print(f"phase_factor = {phase_factor}")
    
    #k = 2 * 1.496e8 * math.pow(10, observer.get_planmag_sun(band) / 5)
    k = 1329
    abs_mag = 5 * math.log10(k / ((body.diam) * math.sqrt(body.albedo)))
    #print(f"abs_mag = {abs_mag}")

    #calculate apparent mag at phase angle and distance from observer given band; this is PLANETARY MAG
    app_mag = abs_mag + 2.5 * math.log10((d_sb**2 * d_bo**2) / (phase_factor * d_os**4))
    #print(f"app_mag of {body} in {band} = {app_mag}")
    
    return app_mag

#returns t/f body is detectable by observers given bands
#body and sun should be instances of Body class
#observers should be a dictionary with keys names and values an instance of the Observer class
#active_bands should be a dictionary with keys observer names and values an array of band names
#mode options are 'tf' or 'info'
#UPDATE: need to address LOS and angular res issues with background
def is_detectable(body, sun, observers, mode='tf', active_bands=None):
    if mode == 'tf':
        for obs in observers:
            if not position_checks(observers[obs], body, sun):
                if active_bands is not None:
                    bands = active_bands[obs]
                else:
                    bands = observers[obs].bands
                for band in bands:
                    if get_app_mag(observers[obs], band, sun, body) <= observers[obs].bands[band][1]:
                        return True
        return False
    
    #return more detailed info: whether or not detected, but also app_mag and band detected in
    elif mode == 'info':
        detected = False
        info = {}
        for obs in observers:
            if not position_checks(observers[obs], body, sun):
                if active_bands is not None:
                    bands = active_bands[obs]
                else:
                    bands = observers[obs].bands
                for band in bands:
                    app_mag = get_app_mag(observers[obs], band, sun, body)
                    if app_mag <= observers[obs].bands[band][1]:
                        #store detection info
                        if obs not in info:
                            info[obs] = []
                        info[obs].append(band)
                        detected = True

        if detected:
            return True, info
        
        return False, None
    
#determine is body is eclipsing sun
def eclipsing_sun(observer, body, sun):
    # #check radii to see if body is closer to observer than sun
    d_bo = dist(observer, body)
    d_os = dist(observer, sun)
    if d_os < d_bo:
        return False
    
    #check if body is eclipsing sun
    theta1 = math.atan(0.5 * sun.diam * 6.6845871226706e-9 / d_os)
    d_sb = dist(body, sun)
    area = 0.25 * math.sqrt((d_bo + d_os + d_sb) * (-d_bo + d_os + d_sb) * (d_bo - d_os + d_sb) * (d_bo + d_os - d_sb)) #heron's
    h = 2 * area / d_os
    theta2 = math.asin(h / d_bo)
    if theta2 <= theta1:
        #print(f"d_ob={d_ob}, d_os={d_os}, d_bs={d_bs}, theta1={theta1}, area={area}, h={h}, theta2={theta2} \n")
        return True
    
    return False

#is the body literally inside the sun
def in_sun(sun, body):
    if dist(sun, body) <= 0.5 * sun.diam * 6.6845871226706e-9:
        return True
    
    return False

#is the phase angle to extreme
def extreme_phase(observer, body, sun):
    d_sb = dist(sun, body)
    d_bo = dist(body, observer)


    max_phase = 160 #HCD
    if math.acos((d_bo**2 + d_sb**2 - dist(observer, sun)**2) / (2 * d_bo * d_sb)) > max_phase * 2 * math.pi / 360 :
        return True
    
    return False

#combine smaller checks into one wrapper
def position_checks(observer, body, sun):
    # if eclipsing_sun(observer, body, sun) or in_sun(sun, body) or extreme_phase(observer, body, sun):
    #     return True
    if in_sun(sun, body) or extreme_phase(observer, body, sun):
        return True
    return False
        
#i hate matplotlib, truly
#adjust coordinates from AU to whatever useless scale pcolormesh sets (as far as I know, you can't change it while preserving heat map)
def plt_sucks(coord, width, pixels):
    return pixels * coord / (2 * width) + pixels / 2
    

    



