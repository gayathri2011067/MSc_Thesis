cd run_files
start=`date +%s`
gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90
./run
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
python3 ../codes/plot.py --trial 14
cd ..
# python3 codes/compare_plots.py 1 "label1" 2 "label2"
#3 "label3" 4 "label4" 5 "label 5" 6 "label 6"
#1 "Model A" 2 "Model B" 3 "Model C" 4 "Model D" 5 "Model E" 6 "Model F" 7 "Model G" 8 "Model H"
