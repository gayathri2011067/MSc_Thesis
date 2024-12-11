cd run_files
start=`date +%s`
gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90
./run
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
python3 ../codes/plot.py --trial 1
cd ..
# python3 codes/compare_plots.py 1 "label1" 2 "label2"
#R_u = 0, R_k = 0
#python3 codes/compare_plots.py 21 "R_u = 1, R_k = 0" 22 "R_u = 0, R_k = 0.3" 23 "R_u = 1, R_k = 0.3" 24 "R_u = 0, R_k = 0"

# python3 codes/compare_plots.py 56 "\alpha_k + \alpha_m \mathrm{from} R_k" 57 "\alpha_k \mathrm{only}" 58 "\alpha_m \mathrm{from} F^{NV}, \xi = 1" 59 "\alpha_m \mathrm{from} F^{NV}, \xi = 0.1"

# python3 codes/compare_plots.py 60 "60" 61 "61" 62 "62" 63 "63" 64 "64" 65 "65"
# python3 codes/compare_plots.py 56 "A" 57 "B" 58 "C" 59 "D"
# python3 codes/compare_plots.py 80 "A: \alpha_k, F^{D}, F^{A} \neq 0;  F^{NV} = 0" 81 "B: \alpha_k \neq 0; \alpha_m = 0" 82 "C: \alpha_k, F^{D} = 0; F^{A}, F^{NV} \neq 0" 84 "E: \alpha_k, F^{D}, F^{A}, F^{NV} \neq 0"
# python3 codes/compare_plots.py 80 "Model A" 81 "Model B" 82 "Model C" 84 "Model E"


# python3 codes/compare_plots.py 84 "\rho \neq f(z); u = f(z)" 85 "\rho, u = f(z)" 
# python3 codes/compare_plots.py 80 "A: \alpha_k, F^{D}, F^{A} \neq 0;  F^{NV} = 0" 81 "B: \alpha_k \neq 0; \alpha_m = 0" 82 "C: \alpha_k, F^{D} = 0; F^{A}, F^{NV} \neq 0" 84 "E: \alpha_k, F^{D}, F^{A}, F^{NV} \neq 0" 85 "\rho, u = f(z)" 


# python3 codes/compare_plots.py 80 "A: \alpha_k, F^{D}, F^{A} \neq 0;  F^{NV} = 0" 82 "C: \alpha_k, F^{D} = 0; F^{A}, F^{NV} \neq 0" 84 "E: \alpha_k, F^{D}, F^{A}, F^{NV} \neq 0" 85 "F: \rho, u = f(z)" 
# python3 codes/compare_plots.py 80 "Conventional model" 82 "Only Vishniac flux" 84 "Vishniac flux + kinetic alpha" 85 "Density stratification"