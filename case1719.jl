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

case_19_xcal_updated_file = CSV.File("mesh_study/Case_19_updated.csv")
case_19_xcal_updated_x = case_19_xcal_updated_file.var"Points:1"
case_19_xcal_updated_h = case_19_xcal_updated_file.var"h"*1e6
case_19_xcal_updated_z_w = case_19_xcal_updated_x/0.1016

case_19_xcal_updated_new_file = CSV.File("mesh_study/Case_19_updated-new.csv")
case_19_xcal_updated_new_x = case_19_xcal_updated_new_file.var"Points:1"
case_19_xcal_updated_new_h = case_19_xcal_updated_new_file.var"h"*1e6
case_19_xcal_updated_new_z_w = case_19_xcal_updated_new_x/0.1016

case_19_xcal_updated_time_stepping_file = CSV.File("mesh_study/Case_19_updated_time_stepping.csv")
case_19_xcal_updated_time_stepping_x = case_19_xcal_updated_time_stepping_file.var"Points:1"
case_19_xcal_updated_time_stepping_h = case_19_xcal_updated_time_stepping_file.var"h"*1e6
case_19_xcal_updated_time_stepping_z_w = case_19_xcal_updated_time_stepping_x/0.1016

case_19_xcal_updated_changed_update_order_file = CSV.File("mesh_study/Case_19_updated_changed_update_order.csv")
case_19_xcal_updated_changed_update_order_x = case_19_xcal_updated_changed_update_order_file.var"Points:1"
case_19_xcal_updated_changed_update_order_h = case_19_xcal_updated_changed_update_order_file.var"h"*1e6
case_19_xcal_updated_changed_update_order_z_w = case_19_xcal_updated_changed_update_order_x/0.1016

case_19_xcal_changed_inlet_to_area_file = CSV.File("mesh_study/Case_19_changed_inlet_to_area.csv")
case_19_xcal_changed_inlet_to_area_x = case_19_xcal_changed_inlet_to_area_file.var"Points:1"
case_19_xcal_changed_inlet_to_area_h = case_19_xcal_changed_inlet_to_area_file.var"h"*1e6
case_19_xcal_changed_inlet_to_area_z_w = case_19_xcal_changed_inlet_to_area_x/0.1016

case_19_xcal_corr_measure_update_ord_file = CSV.File("mesh_study/update_order_changed_correct_measurement_location.csv")
case_19_xcal_corr_measure_update_ord_x = case_19_xcal_corr_measure_update_ord_file.var"Points:1"
case_19_xcal_corr_measure_update_ord_h = case_19_xcal_corr_measure_update_ord_file.var"h"*1e6
case_19_xcal_corr_measure_update_ord_z_w = case_19_xcal_corr_measure_update_ord_x/0.1016

plot(case_17_cfd_z_w, case_17_cfd_h, label = "lan et al. VoF CFD", legend=:bottom)
plot!(case_17_experimental_z_w, case_17_experimental_h, label = "lan et al. Experimental")
plot!(case_17_xcal_z_w, case_17_xcal_h, label = "XCALibre results")
xlabel!("z/W")
ylabel!("h (μm)")

plot(case_19_cfd_z_w, case_19_cfd_h, label = "lan et al. VoF CFD", legend=:bottom, linestyle=:dash)
plot!(case_19_experimental_z_w, case_19_experimental_h, label = "lan et al. Experimental", markershape=:X)
#plot!(case_19_xcal_z_w, case_19_xcal_h, label = "XCALibre results")
plot!(case_19_xcal_updated_z_w, case_19_xcal_updated_h, label = "XCALibre results - old")
plot!(case_19_xcal_updated_new_z_w, case_19_xcal_updated_new_h, label = "XCALibre results - new", linestyle=:dashdot)

plot!(case_19_xcal_updated_time_stepping_z_w, case_19_xcal_updated_time_stepping_h, label = "XCALibre results - fixed time stepping", linestyle=:dashdot)

#plot!(case_19_xcal_updated_changed_update_order_z_w, case_19_xcal_updated_changed_update_order_h, label = "XCALibre results - changed update order", linestyle=:dashdot)
plot!(case_19_xcal_changed_inlet_to_area_z_w, case_19_xcal_changed_inlet_to_area_h, label = "XCALibre results - changed inlet style", linestyle=:dashdot)
plot!(case_19_xcal_corr_measure_update_ord_z_w, case_19_xcal_corr_measure_update_ord_h, label = "XCALibre results - changed update order", linestyle=:dashdot)
xlabel!("z/W")
ylabel!("h (μm)")


plot(case_19_cfd_z_w, case_19_cfd_h, label = "lan et al. VoF CFD", legend=:bottom, linestyle=:dash)
plot!(case_19_xcal_updated_new_z_w, case_19_xcal_updated_new_h, label = "XCALibre results - previous", linestyle=:dashdot)
plot!(case_19_xcal_changed_inlet_to_area_z_w, case_19_xcal_changed_inlet_to_area_h, label = "XCALibre results - changed inlet style")
plot!(case_19_xcal_corr_measure_update_ord_z_w, case_19_xcal_corr_measure_update_ord_h, label = "XCALibre results - changed update order")
xlabel!("z/W")
ylabel!("h (μm)")