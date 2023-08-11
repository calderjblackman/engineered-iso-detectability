#Calder Blackman
#created 7/26/23

#TO-DO tags denote incomplete routines
#UPDATE tags denote code that probably needs to be updated

#math
import math

#astro
from astroquery.jplhorizons import Horizons

#representation of an observatory
#key booleans are:
#i) is_generated: does it actually exist (or, rather, can we query JPL and get a position) 
#ii) on_earth: is it on earth
#iii) is_orbiting: (assumes not on_earth) does it move in space? or is fixed at a grid point
class Observer:
    #bands should be a dictionary with each key the name of the band and each value an array that contains in index 
    #0 the angular res of the band (in degrees)
    #1 the sensitivity of the telescope in the band
    #absmag_sun should be a dictionary with each key the name of the band and each value the absolute magnitude of the sun in the band
    #UPDATE: fov_table should be astropy Table with columns RA_min, RA_max, DEC_min, DEC_max and rows steps
    #required params are is_generated, bands, mag_sun
    #currently expects user to pass correct set of parameters for a given observer type, i.e., not idiot proof
    def __init__(self, is_generated, 
                 #band info
                 bands=None, absmag_sun=None, 
                 #fov constraints
                 fov_table=None, ra_min=0, ra_max=360, dec_min=-90, dec_max=90,
                 #if not generated
                 name=None, origin='@0', epochs=None,
                 #if generated
                 id=None, init_x=None, init_y=None, init_z=None, on_earth=None,
                 #if is_orbiting
                 is_orbiting=False, vect_table=None, a=None, e=None, arg_per=None, long_asc=None, incl=None, 
                 mean_anom=None,  
                 #if write outputs to text file
                 file=None):
        
        #boolean; tracks whether observer is generated, i.e., does this observatory exist and should we query JPL horizons
        self.is_generated = is_generated

        #boolean; tracks if on earth
        self.on_earth = on_earth

        #boolean; tracks if orbiting (if not observer is just fixed at point in cartesian grid)
        self.is_orbiting = is_orbiting

        #store position
        self.x = None
        self.y = None
        self.z = None
    
        #store band info
        self.bands = bands
        self.absmag_sun = absmag_sun

        #calc and store planetary mag of sun in bands
        self.planmag_sun = {}
        for band in self.absmag_sun:
            #calculates planetary magnitude of sun (at 1 AU) in given band
            self.planmag_sun[band] = absmag_sun[band] + 5 * math.log10(1 / 2062650)

        #store fov info
        self.fov_table = fov_table
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max

        #file output
        self.f = None
        if file is not None:
            self.f = file

        #for generated observers
        if self.is_generated:
            #if observer has ID store it
            if id is not None:
                self.id = id

            #initial cartesian coords; with respect to solar system barycenter
            self.x = init_x
            self.y = init_y
            self.z = init_z

            if self.on_earth:
                if self.f is not None:
                    self.f.write(f"OBSERVER created: is_generated=True, on_earth=True, ID={self.id} \n")

            #vector table for orbiting observers
            elif self.is_orbiting:
                #try to load vector table
                self.vect_table = vect_table

                #else do orbital calcs internally based on input orbtal elements
                if self.vect_table is None:
                    #check if necessary orbital params given to calc orbits
                    if not all(param is not None for param in [a, e, arg_per, long_asc, incl, mean_anom]):
                        if self.f is not None:
                            self.f.write(f'ERROR: cannot create Body; 6 orbital params required and <6 given \n')
                        else:  
                            print(f'ERROR: cannot create Body; 6 orbital params required and <6 given')

                    #store orbital elements
                    self.a = a
                    self.e = e
                    self.arg_per = arg_per
                    self.long_asc = long_asc
                    self.incl = incl
                    self.mean_anom = mean_anom

                if self.f is not None:
                    self.f.write(f"OBSERVER created: is_generated=True, is_orbiting=True, ID={self.id} \n")

            #if not is_orbiting and is not on_earth
            else:
                if self.f is not None:
                    self.f.write(f"OBSERVER created: is_generated=True, is_orbiting=False, ID={self.id} \n")
                
        #for existing observatories
        else:
            #check if given and store name for known observers
            if name is None:
                if self.f is not None:
                    self.f.write(f'ERROR: cannot create Observer; is_generated=False but name=None \n')
                else:  
                    print(f'ERROR: cannot create Observer; is_generated=False but name=None')
            self.name = name

            #load vector data from JPL Horizons
            #restricted to name,x,y,z
            #UPDATE?
            if self.on_earth:
                self.vect_horizons = Horizons(id='399', location=origin, epochs=epochs)
                self.vect_table = self.vect_horizons.vectors()['targetname','x','y','z','datetime_str']
            else:
                self.vect_horizons = Horizons(id=self.name, location=origin, epochs=epochs)
                self.vect_table = self.vect_horizons.vectors()['targetname','x','y','z','datetime_str']

            #initial cartesian coords
            self.x = self.vect_table[0]['x']
            self.y = self.vect_table[0]['y']
            self.z = self.vect_table[0]['z']

            if self.f is not None:
                    self.f.write(f"OBSERVER created: is_generated=False, ID={self.name} \n")

        #store initial position for reset after each run
        self.init_x = self.x
        self.init_y = self.y
        self.init_z = self.z

        #display basic info stored about observatory
        if self.f is not None:
            self.f.write(f"    bands: {self.bands.keys()} \n")
            #self.f.write(f"    planmag_sun: {self.planmag_sun} ")

            if self.fov_table is not None:
                self.f.write('    FOV will be determined by fov_table \n')
            elif self.on_earth:
                self.f.write(f"    FOV: ra_min={self.ra_min}, ra_max={self.ra_max}, dec_min={self.dec_min}, dec_max={self.dec_max} \n")
            else:
                self.f.write(f"    FOV is variable, this is a space-based observatory (aka I haven't decided how to handle this) \n")

    #__str__ method; return name if observatory is real (according to JPL), else return ID if exists
    def __str__(self):
        if self.is_generated: 
            if self.id is not None:
                return f"{self.id}"
            else:
                return 'Observer not named: is_generated=True and id=None'

        #for existing observatories
        else: 
            return f"{self.name}"
        
    # ----------------------------------------------------------------------------------- #
    # ------------------------------- Functional Methods  ------------------------------- #
    # ----------------------------------------------------------------------------------- #
    
    #get name/ID if exists
    def get_name(self):
        if self.is_generated: 
            if self.id is not None:
                return f"{self.id}"
            else:
                return 'Observer has no identifier: is_generated=True and id=None'
            
        #for existing observatories    
        else: 
            return f"{self.name}"


    #return triple with x,y,z coordinates
    def get_xyz(self):
        return self.x, self.y, self.z

    #update x,y,z coordinates
    def update_xyz(self, step, iter):
        if self.is_generated:
            if self.on_earth:
                if self.f is not None:
                        self.f.write('patience, this feature is in development \n')
                else:
                    print('patience, this feature is in development')

                    # ---------- TO-DO: need to calc orbits given step ---------- #

            elif self.is_orbiting:
                if self.vect_table is None:
                    if self.f is not None:
                        self.f.write('patience, this feature is in development \n')
                    else:
                        print('patience, this feature is in development')

                    # ---------- TO-DO: need to calc orbits given step ---------- #

                #if self.vect_table is not None
                else:
                    self.x = self.vect_table[iter]['x']
                    self.y = self.vect_table[iter]['y']
                    self.z = self.vect_table[iter]['z']
            
        #if not is_generated          
        else:
            self.x = self.vect_table[iter]['x']
            self.y = self.vect_table[iter]['y']
            self.z = self.vect_table[iter]['z']

    # #check if body is in FOV of telescope at step (which corresponds to row in ephemeris)
    # #UPDATE
    # def in_fov(self, body, step):
    #     #UPDATE: if fov_table was given; expects astropy Table
    #     if self.fov_table is not None:
    #         ra_min = self.fov_table[step]['RA_min']
    #         ra_max = self.fov_table[step]['RA_max']
    #         dec_min = self.fov_table[step]['DEC_min']
    #         dec_max = self.fov_table[step]['DEC_max']
        
    #     ra, dec = body.get_radec(step)

    #     #perform check
    #     if dec >= dec_min and dec <= dec_max and ra >= ra_min and ra <= ra_max:
    #         return True #if in FOV
    #     return False #if not in FOV
    
    # #required param is background, the list of background bodies stored in SolarSystem instance
    # #UPDATE
    # def can_resolve(self, background_vals, body, step):
    #     ra, dec = body.get_radec(step)

    #     for mb in background_vals:
    #         #can't get ra,dec for observer on earth
    #         if not mb.is_earth:
    #             mb_ra, mb_dec = mb.get_radec(step)

    #             approx_dist = math.sqrt((ra - mb_ra)**2 + (dec - mb_dec)**2)

    #             for band in self.bands: 
    #                 if approx_dist <= self.bands[band][0]:
    #                     return False, mb.get_name()
                        
    #     return True, None
    
    #return planmag of sun in given band
    def get_planmag_sun(self, band):
        return self.planmag_sun[band]
    
    #used to make sure type observer is passed
    #lazy, will update later (after the immediate gaping holes in the code are addressed)
    def test():
        pass

    #reset
    def reset(self):
        self.x = self.init_x
        self.y = self.init_y
        self.z = self.init_z
        
    # ----------------------------------------------------------------------------------- #
    # ------------------------------- Methods for Output -------------------------------- #
    # ----------------------------------------------------------------------------------- #

    #set color
    def set_color(self, color): 
        self.color = color
    
    #set marker
    def set_marker(self, marker): 
        self.marker = marker
