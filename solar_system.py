#Calder Blackman
#created 7/28/23

#TO-DO tags denote incomplete routines
#UPDATE tags denote code that probably needs to be updated

#HCD tags denote hardcoded values/parameters

#current limitations to address:
#if you are looking at existing bodies you must pass an observer when they are created in order to get ephem data from JPL
#LOS issues
#angular res issues

#astro
from astroquery.jplhorizons import Horizons

#plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#from simulation
from body import Body
from observer import Observer
from helpers import *

class SolarSystem:
    #background should be dictionary with keys names and values diameter in km
    def __init__(self, name=None, background=None, 
            #output controls
            main_out=None,
            #time controls
            epochs=None, start_time=None, fin_time=None, step=None, step_unit='d'):
        
        #name of this SolarSystem instance
        self.name = name

        #variables used for output control in run()
        self.output_type = None #expects string; controls type of output
        self.anim_type = None #boolean; output animation
        self.fig = None #figure for animation
        self.ax = None #axis for animation
        self.slices_out  =None #boolean; output slices
        self.max_slices = None  #integer; max slices output
        self.slice_plane = None #plane to take slice in
        self.raw_out = None #output raw data
        self.max_raw = None #max number raw files
        self.grid_2d = None #2d grid
        self.grid_3d = None #3d grid
        self.grid_type = None #grid type; None if no grid wanted
        self.app_mag_observer = None #observer for grid
        self.app_mag_band = None #band to do calcs for grid in
        self.grid_width = None #grid width
        self.grid_height = None #grid height
        self.grid_res = None #grid resolution
        self.gridp_diam = None  #diam for grid points
        self.gridp_albedo = None #albedo for grid points
        self.run_name = None #track run name
        self.curr_run_id = 1 #for generating run id's

        #track number of gridpoints
        self.grid_count = 0

         # -------------------------- MAIN OUTPUT FILE SETUP ---------------------------- #

        #name of file that heartbeats and errors are written to
        if main_out is not None:
            self.main_out = main_out
            self.f = open(main_out, 'w')
            self.f.write('\n WHERE ARE THEY HIDING? \n') #WE GON FIND DA ALIENS
            self.f.write('\n created by calder blackman \n')
            self.f.write('\n BSRC internship 2023X \n')
            self.f.write('\n')
        else:
            self.f = None

        # -------------------------------- TIME CONTROLS -------------------------------- #

        #timestep size
        self.step = step

        #set unit for step; deafult is 'd' for days; also 'w','y', etc.; see JPL Horizons docs for full details
        self.step_unit = step_unit

        #store and/or write out epochs to query JPL Horizons with
        if epochs is not None:
            self.epochs = epochs
        else:
            if start_time is None or fin_time is None:
                if self.f is not None:
                    self.f.write('ERROR: incomplete time info; either start_time=None or fin_time=None \n')
                else:
                    print('ERROR: incomplete time info; either start_time=None or fin_time=None')

            self.start_time = start_time
            self.fin_time = fin_time
            #UPDATE (the whole system): default step is in days
            self.epochs = {'start':self.start_time, 'stop':self.fin_time, 'step':f"{self.step}{self.step_unit}"}

        #UPDATE: set max_iter
        this_is_lazy = Horizons(id='A/2017 U1', location='@0', epochs = self.epochs)
        self.max_iter = len(this_is_lazy.vectors()) - 1

        # -------------------------------- GENERAL SETUP -------------------------------- #

        #stores observers with 
        #keys: name/id of observer
        #values: the corresponding Observer instance
        self.observers = {}

        #subset of self.observers
        self.active_observers = None

        #list of active bands; see run() for usage and limitations
        self.active_bands = None

        #expects dictionary with
        #keys: observer name
        #vals: array of with items of type Body
        self.bodies = {}

        #used to generate id's for objects not passed a name
        self.gen_obs_id = 1
        self.gen_body_id = 1

        #setup background bodies; default includes earth; sun always included
        self.background = {}
        if background is not None:
            for name in background:
                if name == '399' or 'earth':
                    self.background[name] = Body(is_generated=False, diam=12756, name=name, epochs=self.epochs, is_earth=True, file=self.f)
                else:
                    self.background[name] = Body(is_generated=False, diam=background[name], name=name, epochs=self.epochs, file=self.f)
        else:
            #add earth as default
            self.background['earth'] = Body(is_generated=False, diam=12756, name='399', epochs=self.epochs, is_earth=True, file=self.f)
        
        #add sun regardless
        self.background['sun'] = Body(is_generated=False, diam=1392000, name='10', epochs=self.epochs, file=self.f)
        
        if self.f is not None:
            self.f.write(f"SOLARSYSTEM created: name={self.name} \n")
           
    def __str__(self):
        if self.name is not None:
            return(self.name)
        return 'SOLARSYSTEM has no identifier'

    #add an Observer instance
    def add_observer(self, is_generated, 
                 #band info
                 bands=None, absmag_sun=None, 
                 #fov constraints
                 fov_table=None, ra_min=0, ra_max=360, dec_min=-90, dec_max=90,
                 #if not generated
                 name=None, origin=None,
                 #if generated
                 id=None, init_x=None, init_y=None, init_z=None, on_earth=None,
                 #if is_orbiting
                 is_orbiting=False, vect_table=None, a=None, e=None, arg_per=None, long_asc=None, incl=None, 
                 mean_anom=None):
        
        #setup identifier
        temp = None
        if name is not None:
            temp = str(name)
        elif id is not None:
            temp = str(id)
        else:
            temp = 'obs'+str(self.gen_obs_id)
            self.gen_obs_id += 1

        #if, elif series prevents observers being duplicated
        if type(temp) is str:
            if temp not in self.observers:
                #instantiate object
                self.observers[temp] = Observer(is_generated=is_generated,
                        #band info
                        bands=bands, absmag_sun=absmag_sun, 
                        #fov constraints
                        fov_table=fov_table, ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max,
                        #if not generated
                        name=name, origin=origin, epochs=self.epochs,
                        #if generated
                        id=temp, init_x=init_x, init_y=init_y, init_z=init_z, on_earth=on_earth,
                        #if is_orbiting
                        is_orbiting=is_orbiting, vect_table=vect_table, a=a, e=e, arg_per=arg_per, long_asc=long_asc, incl=incl, 
                        mean_anom=mean_anom,  
                        #if write outputs to text file
                        file=self.f)
                
                #store in SolarSystem
                #self.active_observers[temp] = self.observers[temp]
                self.bodies[temp] = []

        elif type(temp) is int:
            if int >= self.gen_obs_id:
                #instantiate object
                self.observers[temp] = Observer(is_generated=is_generated,
                        #band info
                        bands=bands, absmag_sun=absmag_sun, 
                        #fov constraints
                        fov_table=fov_table, ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max,
                        #if not generated
                        name=name, origin=origin, epochs=self.epochs,
                        #if generated
                        id=temp, init_x=init_x, init_y=init_y, init_z=init_z, on_earth=on_earth,
                        #if is_orbiting
                        is_orbiting=is_orbiting, vect_table=vect_table, a=a, e=e, arg_per=arg_per, long_asc=long_asc, incl=incl, 
                        mean_anom=mean_anom,  
                        #if write outputs to text file
                        file=self.f)
                
                #store in SolarSystem
                #self.active_observers[temp] = self.observers[temp]
                self.bodies[temp] = []

    #add a Body instance     
    def add_body(self, is_generated,
                 #physical characteristics
                 diam=None, albedo=None, 
                 #observer, should be Observer instance
                 observer=None,
                 #if not generated
                 name=None, origin=None, is_earth=False,
                 #if generated
                 id=None, init_x=None, init_y=None, init_z=None, 
                 #if is_orbiting
                 is_orbiting=False, vect_table=None, ephem_table=None, a=None, e=None, arg_per=None, long_asc=None, incl=None, 
                 mean_anom=None):
        
        #handle when observer is none; all bodies store with key 'NONE
        if observer is None:
            temp_observer = 'NONE'
            
            if self.f is not None:
                self.f.write('WARNING: observer=None in add_body() call \n')
            else:
                print('WARNING: observer=None in add_body() call')
        else:
            temp_observer = observer.get_name()

        #setup identifier
        temp = None
        if name is not None:
            temp = str(name)
        elif id is not None:
            temp = str(id)
        else:
            temp = 'body'+str(self.gen_body_id)
            self.gen_body_id += 1

        #instantiate object
        self.bodies[temp_observer].append(Body(is_generated=is_generated,
                 #physical characteristics
                 diam=diam, albedo=albedo, 
                 #observer, should be Observer instance
                 observer=observer,
                 #if not generated
                 name=name, origin=origin, epochs=self.epochs, is_earth=is_earth,
                 #if generated
                 id=temp, init_x=init_x, init_y=init_y, init_z=init_z, 
                 #if is_orbiting
                 is_orbiting=is_orbiting, vect_table=vect_table, ephem_table=ephem_table, a=a, e=e, arg_per=arg_per, long_asc=long_asc, 
                 incl=incl, mean_anom=mean_anom, 
                 #if write outputs to text file
                 file=self.f))
        
    #reset observers
    def clear_observers(self):
        self.observers = {}
        self.active_observers = None
        self.active_bands = None
        self.gen_obs_id = 1

    # ----------------------------------------------------------------------------------- #
    # ---------------------------------- Run() Method ----------------------------------- #
    # ----------------------------------------------------------------------------------- #

    #call this function to run; allows various runs and outputs with same SolarSystem instance
    #expects init_time, fin_time in form 'YEAR-MONTH-DAY'; for example:'2016-01-01'
    #expects step as integer
    def run(self, run_name=None, output_type=None,
            #output controls
            anim_type=None, slices_out=False, max_slices=5, slice_plane=None, raw_out=False, max_raw=5,
            #grid controls
            grid_type=None, grid_width=None, grid_height=None, grid_res=None, 
            gridp_diam=None, gridp_albedo=None,
            #for output_type=='app_mag'
            app_mag_observer=None, app_mag_band=None,
            #for output_type=='observability'; 
            #active_observers should be array with names of observers
            #active_bands should be a dictionary with keys observer names and values an array of band names
            active_observers=None, active_bands=None):
        
        if self.f is not None:
            self.f.write('run(): ROUTINE BEGINNING \n')

        #reset everything before new run
        self.reset()

        # ------------------------------- OUTPUT CONTROLS ------------------------------- #

        #string; controls type of outputs
        #options (for now) are:
        #'observability': this gives an output showing t/f whether bodies are detected by active observers in active bands
        #'app_mag': designed to work with grid of bodies, not individual bodies; this gives a heat map of apparent magnitudes calculated
        if output_type is None:
            if self.f is not None:
                self.f.write('ERROR: no output type specified; output_type=None')
            else:
                print('ERROR: no output type specified; output_type=None')
        self.output_type = output_type

        #controls if animation is output; expects string, either '2d' or '3d'
        self.anim_type = anim_type 
            
        #boolean; output slices (equivalent to frames in animation); these are 2D slices on z=0 plane (for now)
        self.slices_out = slices_out

        #integer; max number of slices output; default is 5
        self.max_slices = max_slices

        #plane to take slice in, in terms of height (no fancy planes, yet)
        self.slice_plane = slice_plane

        #boolean; output raw data arrays
        self.raw_out = raw_out

        #nmax number of raw files to output
        self.max_raw = max_raw

        #controls if grid is used; expects string, either '2d' or '3d'
        self.grid_type = grid_type

        #observer to use for grid calcs; default is earth
        #if not default must be type Observer
        self.app_mag_observer = app_mag_observer

        #band to use for grid calcs
        #expects string; only required if output type is 'app_mags'
        #should be a band of the app_mag_observer, obviously
        self.app_mag_band = app_mag_band

        #width of grid in both x and y directions in AU
        #this is radius, i.e., distance either side of sun to build grid
        self.grid_width = grid_width

        #height of grid in AU
        #this is radius, i.e., distance either side of sun to build grid
        self.grid_height = grid_height

        #resolution of grid; dist between bodies in KILOMETERS
        self.grid_res = grid_res

        #convert res to AU
        self.grid_res = self.grid_res * 6.6845871226706e-9

        #physical parameters for grid points
        self.gridp_diam = gridp_diam
        self.gridp_albedo = gridp_albedo

        #give the run a name
        if run_name is None:
            self.run_name = self.curr_run_id
            self.curr_run_id += 1
        else:
            self.run_name = run_name

        # --------------------------------- BUILD GRID ---------------------------------- #

        if self.grid_type is not None:
            #check necessary parameters are defined
            if self.gridp_diam is None or self.gridp_albedo is None:
                if self.f is not None:
                    self.f.write(f"ERROR: gridp_diam={self.gridp_diam} and/or gridp_albedo={self.gridp_albedo} \n")
                else:
                    print(f"ERROR: gridp_diam={self.gridp_diam} and/or gridp_albedo={self.gridp_albedo}")

            #build 2d grid
            if self.grid_type == '2d':
                self.grid_2d = []
                x = -self.grid_width
                while x < self.grid_width:
                    temp_row = []
                    y = -self.grid_width
                    while y < self.grid_width:
                        temp_row.append(Body(is_generated=True,
                                        #physical characteristics
                                        diam=self.gridp_diam, albedo=self.gridp_albedo, 
                                        #if generated
                                        id=f"BODY({x},{y},{0})", init_x=x, init_y=y, init_z=0,
                                        #file
                                        file=self.f, write=False))
                        self.grid_count += 1
                        y += self.grid_res

                    self.grid_2d.append(temp_row)
                    x += self.grid_res

            #build 3d grid
            elif self.grid_type == '3d':
                self.grid_3d = []
                x = -self.grid_width
                while x < self.grid_width:
                    temp_row = []
                    y = -self.grid_width
                    while y < self.grid_width:
                        temp_col = []
                        z = -self.grid_height
                        while z < self.grid_height:
                            temp_col.append(Body(is_generated=True,
                                            #physical characteristics
                                            diam=self.gridp_diam, albedo=self.gridp_albedo, 
                                            #if generated
                                            id=f"BODY({x},{y},{z})", init_x=x, init_y=y, init_z=z,
                                            #file
                                            file=self.f, write=False))
                            self.grid_count += 1
                            z += self.grid_res
                        temp_row.append(temp_col)
                        y += self.grid_res
                    self.grid_3d.append(temp_row)
                    x += self.grid_res

            else:
                if self.f is not None:
                    self.f.write('ERROR in run(): grid_type paramter must be None, 2d, or 3d \n')
                else:
                    print('ERROR in run(): grid_type parameter must be None, 2d, or 3d')

            if self.f is not None:
                self.f.write(f"run(): GRID constructed; grid_type={self.grid_type}, grid_count={self.grid_count} \n")

            #set grid observer and band
            if app_mag_observer is None or app_mag_band is None:
                if self.f is not None:
                    self.f.write('ERROR: app_mag_observer=None and/or app_mag_band=None \n')
                else:
                    print('ERROR: app_mag_observer=None and/or app_mag_band=None')
       
                
            if type(app_mag_observer) == str:
                if app_mag_observer in self.observers:
                    self.app_mag_observer = self.observers[app_mag_observer]
                else:
                    if self.f is not None:
                        self.f.write('ERROR: app_mag_observer is str but does not exist \n')
                    else:
                        print('ERROR: app_mag_observer is str but does not exist')

            else:
                self.app_mag_observer = app_mag_observer

                try:   
                    self.app_mag_observer.test()
                except:         
                    if self.f is not None:
                        self.f.write('ERROR: if app_mag_observer is not existing it must be new type Observer \n')
                    else:
                        print('ERROR: if app_mag_observer is not existing it must be new type Observer')

            self.app_mag_band = app_mag_band

        # --------------------------------- SET ACTIVES --------------------------------- #
        
        if self.output_type == 'observability':
            #set active observers for run
            if active_observers is not None:
                if len(active_observers) == 0:
                    if self.f is not None:
                        self.f.write('WARNING: no active observers; active_observers=[] \n')
                    else:
                        print('WARNING: no active observers; active_observers=[]')

                self.active_observers = {}
                for obs in active_observers:
                    self.active_observers[obs] = self.observers[obs]

            #set active bands for run; only allowed if a single observer is active (for now)
            if active_bands is not None:
                if len(active_bands) == 0:
                    if self.f is not None:
                        self.f.write('WARNING: no active bans; active_bands=[] \n')
                    else:
                        print('WARNING: no active band; active_bands=[]')

                self.active_bands = active_bands

        # -------------------------- MAIN LOOP (no animation) --------------------------- #
        
        if self.anim_type is None:
            #if writing slices
            if self.slices_out:
                curr_slice = 1
                
                #controls how often a slice is written
                slice_freq = self.max_iter // self.max_slices

            if self.raw_out:
                curr_raw = 1

                #controls how often raw file is written
                raw_freq = self.max_iter // self.max_raw

            # ----- (actual) MAIN LOOP ----- #
            for iter in range(self.max_iter):
                #write out heartbeat
                if iter % 10 == 0 or iter == 0:
                    if self.f is not None:
                        self.f.write(f"HEARTBEAT: iteration={iter} \n")

                # ----- UPDATE ALL POSITIONS ----- #
                self.update(iter=iter)

                # ----- SLICE OUTPUT ---- #
                if self.slices_out:
                    #limit number of slices output (slice_freq should handle, for later use); don't want to write a bajillion files
                    if curr_slice > self.max_slices:
                        curr_slice = 1

                    #only output (roughly) min number of slices requested; don't keep overwriting
                    if iter % slice_freq == 0:
                        #for app_mag output; this is grid of app mags
                        if self.output_type == 'app_mag':
                            if self.grid_type is not None:
                                if self.grid_type == '2d':
                                    app_mags = []

                                    for row in self.grid_2d:
                                        app_mag_row = []

                                        for body in row:
                                            #HCD in position_checks()
                                            if position_checks(self.app_mag_observer, body, self.background['sun']):
                                                app_mag_row.append(0)
                                            else:
                                                #do app mag calcs
                                                app_mag_row.append(get_app_mag(self.app_mag_observer, self.app_mag_band, self.background['sun'], body))

                                        app_mags.append(app_mag_row)

                                    #create heatmap plot
                                    fig, ax = plt.subplots()
                                    hmap = ax.pcolormesh(app_mags)

                                    #plot background bodies
                                    for name in self.background:
                                        if np.linalg.norm(np.array([self.background[name].get_xyz()])) < self.grid_width:
                                            ax.plot(plt_sucks(self.background[name].get_xyz()[1], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                                    plt_sucks(self.background[name].get_xyz()[0], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                                    color='black', marker='o')

                                    #configer and save plot        
                                    ax.set_xticks([])
                                    ax.set_yticks([])
                                    fig.colorbar(hmap)
                                    ax.set_title(f"{self.run_name}_slice{curr_slice}_appmag2d (iter={iter})")  
                                    plt.savefig(f"outputs/{self.run_name}_slice{curr_slice}_appmag2d")

                                elif self.grid_type == '3d':
                                    pass

                                    # ---------- TO-DO: need to do 3d grid get_app_mag output ---------- #

                                    
                            #grid_type = None; no using a grid
                            else:
                                pass

                                # ---------- TO-DO: need to do observability slice output ---------- #


                        elif self.output_type == 'observability':

                            # ---------- TO-DO: need to do observability slice output ---------- #
                            #need consider both grid and not

                            if self.grid_type is not None:
                                if self.grid_type == '2d':
                                    detected = []

                                    if self.active_observers is not None:
                                        for row in self.grid_2d:
                                            detected_row = []

                                            for body in row:
                                                if is_detectable(body, self.background['sun'], self.active_observers, 
                                                                 active_bands=self.active_bands):
                                                    detected_row.append(1)
                                                else:
                                                    detected_row.append(0)
                                            detected.append(detected_row)
                                    
                                    #use all observers, not active_observers
                                    else:
                                        for row in self.grid_2d:
                                            detected_row = []

                                            for body in row:
                                                if is_detectable(body, self.background['sun'], self.observers, 
                                                                 active_bands=self.active_bands):
                                                    detected_row.append(1)
                                                else:
                                                    detected_row.append(0)
                                            detected.append(detected_row)           
                                                    
                                    #create heatmap plot
                                    fig, ax = plt.subplots()
                                    hmap = ax.pcolormesh(detected)

                                    #plot background bodies
                                    for name in self.background:
                                        if np.linalg.norm(np.array([self.background[name].get_xyz()])) < self.grid_width:
                                            ax.plot(plt_sucks(self.background[name].get_xyz()[1], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                                    plt_sucks(self.background[name].get_xyz()[0], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                                    color='black', marker='o')

                                    #configer and save plot        
                                    ax.set_xticks([])
                                    ax.set_yticks([])
                                    fig.colorbar(hmap)
                                    ax.set_title(f"{self.run_name}_slice{curr_slice}_obs2d (iter={iter})")  
                                    plt.savefig(f"outputs/{self.run_name}_slice{curr_slice}_obs2d")

                                elif self.grid_type == '3d':
                                    pass

                            #grid_type = None; no using a grid
                            else:
                                pass

                                # ---------- TO-DO: need to do observability slice output ---------- #

                        curr_slice += 1

                # ----- RAW FILE OUTPUT ----- #
                if self.raw_out:
                    #limit number of slices output (slice_freq should handle, for later use); don't want to write a bajillion files
                    if curr_raw > self.max_raw:
                        curr_raw = 1

                    #only output (roughly) min number of slices requested; don't keep overwriting
                    if iter % raw_freq == 0:
                        #for app_mag output; this is grid of app mags
                        if self.output_type == 'app_mag':
                            if self.grid_type is not None:
                                if self.grid_type == '2d':
                                    pass

                                elif self.grid_type == '3d':
                                    raw = []

                                    #build 3d array of apparent magnitudes
                                    for row in self.grid_3d:
                                        raw_row = []
                                        for col in row:
                                            raw_col = []
                                            for body in col:
                                                #HCD in position_checks()
                                                if position_checks(self.app_mag_observer, body, self.background['sun']):
                                                    app_mag = 0
                                                else:
                                                    #do app mag calcs
                                                    app_mag = get_app_mag(self.app_mag_observer, self.app_mag_band, self.background['sun'], body)
                                                raw_col.append(app_mag)
                                            raw_row.append(raw_col)
                                        raw.append(raw_row)
                                    
                                    #reshape array
                                    raw = np.array(raw)
                                    np.save(f"raw_outputs/{self.run_name}_raw{curr_raw}_appmag3d", raw)

                            #grid_type = None; no using a grid
                            else:
                                pass

                                # ---------- TO-DO: need to do observability slice output ---------- #


                        elif self.output_type == 'observability':
                            if self.grid_type == '2d':
                                pass

                            # ---------- TO-DO: need to do observability slice output ---------- #
                            #need consider both grid and not

                            elif self.grid_type == '3d':
                                pass

                            #grid_type = None; no using a grid
                            else:
                                pass

                                # ---------- TO-DO: need to do observability slice output ---------- #

                        curr_raw += 1

        # -------------------------- MAIN LOOP (w/ animation) --------------------------- #

        else:
            if self.anim_type == '2d':
                self.fig = plt.figure()
                self.ax = plt.axes()
            elif self.anim_type == '3d':
                self.fig = plt.figure()
                self.ax = plt.axes(projection='3d')
            else:
                if self.f is not None:
                    self.f.write('ERROR: invalid anim_type \n')
                else:
                    print('ERROR: invalid anim_type')

            anim = FuncAnimation(self.fig, self.animate, frames=self.max_iter)
            anim.save(f"outputs/{self.run_name}_anim_appmag2d.mp4")

        # -------------------------------- FINISH run() --------------------------------- #

        if self.f is not None:
            self.f.write('run(): ROUTINE COMPLETED')

    # ----------------------------------------------------------------------------------- #
    # ------------------------------- Functional Methods  ------------------------------- #
    # ----------------------------------------------------------------------------------- #
    
    #update all bodies and observers in solar system
    def update(self, iter):
        for obs in self.observers:
            self.observers[obs].update_xyz(step=self.step, iter=iter)

            for body in self.bodies[obs]:
                body.update_xyz(step=self.step, iter=iter)
        
        for mb in self.background:
            self.background[mb].update_xyz(step=self.step, iter=iter)

    #reset all bodies and observers in solar system
    def reset(self):
        for obs in self.observers:
            self.observers[obs].reset()

            for body in self.bodies[obs]:
                body.reset()
        
        for mb in self.background:
            self.background[mb].reset()

        if self.f is not None:
            self.f.write('reset() \n')

    # ----------------------------------------------------------------------------------- #
    # ------------------------------- Methods for Output -------------------------------- #
    # ----------------------------------------------------------------------------------- #

    #function to be called by FuncAnimation
    def animate(self, iter):
        #clear plot
        self.ax.clear()

        #update all bodies
        self.update(iter=iter)

        if self.output_type == 'app_mag':
            if self.grid_type is not None:
                if self.grid_type == '2d':
                    app_mags = []

                    for row in self.grid_2d:
                        app_mag_row = []

                        for body in row:
                            #HCD in position_checks()
                            if position_checks(self.app_mag_observer, body, self.background['sun']):
                                app_mag_row.append(0)
                            else:
                                #do app mag calcs
                                app_mag_row.append(get_app_mag(self.app_mag_observer, self.app_mag_band, self.background['sun'], body))

                        app_mags.append(app_mag_row)

                    #create heatmap plot
                    hmap = self.ax.pcolormesh(app_mags)

                    #plot background bodies
                    for name in self.background:
                        if np.linalg.norm(np.array([self.background[name].get_xyz()])) < self.grid_width:
                            self.ax.plot(plt_sucks(self.background[name].get_xyz()[1], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                    plt_sucks(self.background[name].get_xyz()[0], self.grid_width, 2 * self.grid_width // self.grid_res), 
                                    color='black', marker='o')

                    #configer and save plot        
                    self.ax.set_xticks([])
                    self.ax.set_yticks([])
                    if iter == self.max_iter:
                        self.fig.colorbar(hmap)
                    self.ax.set_title(f"{self.run_name}_anim_appmag2d")  

                elif self.grid_type == '3d':
                    pass

                    # ---------- TO-DO: need to do 3d grid app_mag output ---------- #

                    
            #grid_type = None; no using a grid
            else:
                pass

                # ---------- TO-DO: need to do observability animation output ---------- #


        elif self.output_type == 'observability':
            pass

            # ---------- TO-DO: need to do observability animation output ---------- #
            #need consider both grid and not
    
    #close main output file
    def close_f(self):
        if self.f is not None:
            self.f.write('\n')
            self.f.write('\n CLOSING FILE \n')
            self.f.write('\n ITS TIME TO SEE: \n')
            self.f.write('\n WHERE ARE THEY HIDING? \n')
            self.f.close()
