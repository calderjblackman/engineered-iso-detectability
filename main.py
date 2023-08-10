#sample run of SolarSytem code
#best of luck

from solar_system import SolarSystem

ss = SolarSystem(name='ss1', main_out='ss1_out.txt', start_time='2016-01-01', fin_time='2020-01-01', step=21)

#pan-starrs
pst_bands = {'g':(0.06, 23.3), 'r':(0.06, 23.2), 'i':(0.06, 23.1), 'z':(0.06, 22.3), 'y':(0.06, 21.4)}
pst_absmag_sun = {'g':5.14, 'r':4.53, 'i':4.18, 'z':4.02, 'y':3.99}
ss.add_observer(is_generated=False, bands=pst_bands, absmag_sun=pst_absmag_sun, name='pan_starrs', origin='@0', on_earth=True)

# ss.run(run_name='run1', output_type='app_mag', slices_out=True, max_slices=2, 
#        #grid controls
#        grid_type='2d', grid_observer='pan_starrs', grid_band='g', grid_width=3.5, res=10000000, gridp_diam=1, gridp_albedo=0.6)

ss.run(run_name='run2', output_type='app_mag', anim_type='2d', grid_type='2d', grid_observer='pan_starrs', 
       grid_band='g', grid_width=3.5, res=10000000, gridp_diam=1, gridp_albedo=0.6)

ss.close_f()