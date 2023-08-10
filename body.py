#Calder Blackman
#created 7/28/23

#TO-DO tags denote incomplete routines
#UPDATE tags denote code that probably needs to be updated

#astro
from astroquery.jplhorizons import Horizons

#representation of body/point of interest
#key booleans are:
#i) is_generated: does it actually exist (or, rather, can we query JPL and get a position) 
#ii) is_orbiting: does it move in space? or is fixed at a grid point
#currently expects input vector/ephemeris tables to be in same form as tables from astroquery.jplhorizons
class Body:
    #required params are is_generated
    #currently expects user to pass correct set of parameters for a given body type, i.e., not idiot proof
    def __init__(self, is_generated,
                 #physical characteristics
                 diam=None, albedo=None, 
                 #observer, should be Observer instance
                 observer=None,
                 #if not generated
                 name=None, origin='@0', epochs=None, is_earth=False,
                 #if generated
                 id=None, init_x=None, init_y=None, init_z=None, 
                 #if is_orbiting
                 is_orbiting=False, vect_table=None, ephem_table=None, a=None, e=None, arg_per=None, long_asc=None, incl=None, 
                 mean_anom=None, 
                 #if write outputs to text file
                 file=None, write=True):
        
        #tracks whether body is generated, i.e., does this body exist and should we query JPL horizons
        self.is_generated = is_generated
        
        #physical characteristics
        self.diam = diam  #diam in km
        self.albedo = albedo

        #store observer; expects Observer instance
        self.observer=observer

        #boolean; tracks if object was detected in previous timestep
        self.detected = False

        #store position
        self.x = None
        self.y = None
        self.z = None

        #file output for Body objects
        self.f = None
        if file is not None:
            self.f = file

        #for bodies that are generated
        if self.is_generated:
            #if object has ID store it
            if id is not None:
                self.id = id

            #initial cartesian coords; with respect to solar system barycenter
            self.x = init_x
            self.y = init_y
            self.z = init_z

            #is the object orbiting, i.e., is it a gridpoint or a dynamic body
            self.is_orbiting = is_orbiting
            if self.is_orbiting:
                #try to load vector table
                self.vect_table = vect_table

                #try to load ephemeris table
                self.ephem_table = ephem_table

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

                #won't write unless file passed to Body object when instantiated
                if self.f is not None and write:
                    self.f.write(f"BODY created: is_generated=True, is_orbiting=True, ID={self.id} \n")

            else:
                #won't write unless file passed to Body object when instantiated
                if self.f is not None and write:
                    self.f.write(f"BODY created: is_generated=True, is_orbiting=False, ID={self.id} \n")

        #for bodies that are not generated
        else:
            #check if given and store name for known objects
            if name is None:
                if self.f is not None:
                    self.f.write(f'ERROR: cannot create Body; is_generated=False but name=None \n')
                else:  
                    print(f'ERROR: cannot create Body; is_generated=False but name=None')
            self.name = name

            #kinda self explanatory; is this earth? 
            self.is_earth = is_earth

            #load vector data from JPL Horizons
            #restricted to name,x,y,z
            self.vect_horizons = Horizons(id=f"{self.name}", location=origin, epochs=epochs)
            self.vect_table = self.vect_horizons.vectors()['targetname','x','y','z','datetime_str']

            #initial cartesian coords
            self.x = self.vect_table[0]['x']
            self.y = self.vect_table[0]['y']
            self.z = self.vect_table[0]['z']

            #fyi you can't load ephem data for yourself
            if not self.is_earth and self.observer is not None:
                #load ephemeris data from JPL Horizons
                #quantities=1 loads RA and DEC only; further restricted to name,RA,DEC
                self.ephem_horizons = Horizons(id=f"{self.name}", location=self.observer.get_name(), epochs=epochs)
                self.ephem_table = self.ephem_horizons.ephemerides(quantities=1)['targetname','RA','DEC']

            #won't write unless file passed to Body object when instantiated
            if self.f is not None and write:
                self.f.write(f"BODY created: is_generated=False, name={self.name} \n")

        #store initial position for reset after each run
        self.init_x = self.x
        self.init_y = self.y
        self.init_z = self.z

    #__str__ method; return name if object is real (according to JPL), else return ID if exists
    def __str__(self):
        if self.is_generated: 
            if self.id is not None:
                return f"{self.id}"
            else:
                return 'Body has no identifier: is_generated=True and id=None'
            
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
                return 'Body has no identifier: is_generated=True and id=None'
            
        else: 
            return f"{self.name}"
        
    #return triple with x,y,z coordinates
    def get_xyz(self):
        return self.x, self.y, self.z
    
    #update x,y,z coordinates
    def update_xyz(self, step, iter):
        if self.is_generated:
            if self.is_orbiting:
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

    #return tuple with RA,DEC given step (which corresponds to row in ephemeris table)
    def get_radec(self, step=None):
        if self.is_generated:
            if self.ephem_table is not None:
                #UPDATE: expects astropy Table
                return self.ephem_table[step]['RA'], self.ephem_table[step]['DEC']
            else:
                if self.f is not None:
                        self.f.write('patience, this feature is in development \n')
                else:
                    print('patience, this feature is in development')

                # ---------- TO-DO: need to calc ra,dec given x,y,z---------- #
            
        #if not is_generated
        else:
            return self.ephem_table[step]['RA'], self.ephem_table[step]['DEC']
        
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
