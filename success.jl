using Plots

cd_1_success_case=[8]
cd_1_success_angl=[90]
cd_1_bad_sol_case=[1,6,7]
cd_1_bad_sol_angl=[5,90,90]
cd_1_sol_fail_cas=[2,3,4]
cd_1_sol_fal_angl=[5,5,5]
cd_1_par_fail_cas=[5]
cd_1_par_fal_angl=[90]
cd_1_angle=[5,5,5,5,90,90,90,90]


cd_2_par_fail_cas=[9,10,11,12,13,14,15,16,18]
cd_2_success_case=[17,19]
cd_2_angle=[60, 60, 60, 60,60,60,60,60,30,90,30,90]


scatter![cd_1_bad_sol_angl, cd_1_bad_sol_case]
scatter[cd_1_success_case, cd_1_success_angl]