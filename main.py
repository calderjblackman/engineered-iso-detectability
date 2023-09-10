#sample run of SolarSytem code
#best of luck

from solar_system import SolarSystem

ss = SolarSystem(name='ss1', main_out='ss1_out.txt', start_time='2016-01-01', fin_time='2020-01-01', step=21)
#ss.clear_observers()

#pan-starrs
pst_bands = {'g':(0.06, 23.3), 'r':(0.06, 23.2), 'i':(0.06, 23.1), 'z':(0.06, 22.3), 'y':(0.06, 21.4)}
pst_absmag_sun = {'g':5.14, 'r':4.53, 'i':4.18, 'z':4.02, 'y':3.99}
ss.add_observer(is_generated=False, bands=pst_bands, absmag_sun=pst_absmag_sun, name='pan-starrs', origin='@0', on_earth=True)

ss.run(run_name='run1', output_type='app_mag', slices_out=True, max_slices=2, 
       grid_type='2d', app_mag_observer='pan-starrs', app_mag_band='g', grid_width=3.5, grid_res=10000000, gridp_diam=1, gridp_albedo=0.6)

# ss.run(run_name='run2', output_type='app_mag', anim_type='2d', grid_type='2d', app_mag_observer='pan-starrs', 
#        grid_band='g', grid_width=3.5, grid_res=10000000, gridp_diam=1, gridp_albedo=0.6)

#???
# ss.run(run_name='run3', output_type='observability', slices_out=True, max_slices=2, grid_type='2d', grid_width=3.5, 
#        grid_res=10000000, gridp_diam=1, gridp_albedo=0.6)

jsn_bands = {'v':(0.1, 23.0)}
jsn_absmag_sun = {'v': 4.81}
ss.add_observer(is_generated=False, bands=jsn_bands, absmag_sun=jsn_absmag_sun, name='johnson', origin='@0', on_earth=True)

# ss.run(run_name='run4', output_type='app_mag', raw_out=True, max_raw=5, grid_type='3d', 
#        grid_width=3.5, grid_height=3.5, grid_res=20000000, gridp_diam=1, gridp_albedo=0.06, app_mag_observer='johnson', app_mag_band='v')

ss.close_f()