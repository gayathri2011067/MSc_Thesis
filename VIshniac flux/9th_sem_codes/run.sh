cd run_files
start=`date +%s`
gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90
./run
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
python3 ../codes/plot.py --trial 17

cd ..
# python3 codes/compare_plots.py 1 "label1" 2 "label2"
#R_u = 0, R_k = 0
#python3 codes/compare_plots.py 21 "R_u = 1, R_k = 0" 22 "R_u = 0, R_k = 0.3" 23 "R_u = 1, R_k = 0.3" 24 "R_u = 0, R_k = 0"
#python3 codes/compare_plots.py 10 "u = f(z), \rho \neq f(z)" 11 "\rho = f(z), u \neq f(z)" 12 "u = f(z), \rho = f(z)" 
