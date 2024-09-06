cd run_files
start=`date +%s`
gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90
./run
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
python3 ../codes/plot.py --trial 2
cd ..
# python3 codes/compare_plots.py 1 "label1" 2 "label2"