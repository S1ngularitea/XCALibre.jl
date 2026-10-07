using CSV
using Plots

case_17_file = CSV.File("mesh_study/case 17.csv")
case_17_cfd_x = case_17_file.x[1:32]
case_17_cfd_h = case_17_file.var" y"[1:32]
case_17_cfd_z_w = case_17_cfd_x/101.6
case_17_experimental_x = case_17_file.x[33:53]
case_17_experimental_h = case_17_file.var" y"[33:53]
case_17_experimental_z_w = case_17_experimental_x/101.6

case_17_xcal_file = CSV.File("mesh_study/c17.csv")
case_17_xcal_x = case_17_xcal_file.var"Points:1"
case_17_xcal_h = case_17_xcal_file.var"h"*1e6
case_17_xcal_z_w = case_17_xcal_x/0.1016

case_19_file = CSV.File("mesh_study/case 19.csv")
case_19_cfd_x = case_19_file.x[1:30]
case_19_cfd_h = case_19_file.var" y"[1:30]
case_19_cfd_z_w = case_19_cfd_x/101.6
case_19_experimental_x = case_19_file.x[31:52]
case_19_experimental_h = case_19_file.var" y"[31:52]
case_19_experimental_z_w = case_19_experimental_x/101.6

case_19_xcal_file = CSV.File("mesh_study/c19 cd2_2.csv")
case_19_xcal_x = case_19_xcal_file.var"Points:1"
case_19_xcal_h = case_19_xcal_file.var"h"*1e6
case_19_xcal_z_w = case_19_xcal_x/0.1016

plot(case_17_cfd_z_w, case_17_cfd_h, label = "lan et al. VoF CFD", legend=:bottom)
plot!(case_17_experimental_z_w, case_17_experimental_h, label = "lan et al. Experimental")
plot!(case_17_xcal_z_w, case_17_xcal_h, label = "XCALibre results")
xlabel!("z/W")
ylabel!("h (μm)")

plot(case_19_cfd_z_w, case_19_cfd_h, label = "lan et al. VoF CFD", legend=:bottom)
plot!(case_19_experimental_z_w, case_19_experimental_h, label = "lan et al. Experimental")
plot!(case_19_xcal_z_w, case_19_xcal_h, label = "XCALibre results")
xlabel!("z/W")
ylabel!("h (μm)")