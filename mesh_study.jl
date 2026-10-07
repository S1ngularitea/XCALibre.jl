using CSV
using Plots
## Rread Data
CD1_1_File = CSV.File("mesh_study/c8 cd1_1.csv")
CD1_1_h = CD1_1_File.h * 1e6
CD1_1_z = CD1_1_File.var"Points:1"
CD1_1_z_w = CD1_1_z/0.61

CD1_2_File = CSV.File("mesh_study/c8 cd1_2.csv")
CD1_2_h = CD1_2_File.h * 1e6
CD1_2_z = CD1_2_File.var"Points:1"
CD1_2_z_w = CD1_2_z/0.61

CD1_3_File = CSV.File("mesh_study/c8 cd1_3.csv")
CD1_3_h = CD1_3_File.h * 1e6
CD1_3_z = CD1_3_File.var"Points:1"
CD1_3_z_w = CD1_3_z/0.61


CD2_1_File = CSV.File("mesh_study/c19 cd2_1.csv")
CD2_1_h = CD2_1_File.h * 1e6
CD2_1_z = CD2_1_File.var"Points:1"
CD2_1_z_w = CD2_1_z/0.61

CD2_2_File = CSV.File("mesh_study/c19 cd2_2.csv")
CD2_2_h = CD2_2_File.h * 1e6
CD2_2_z = CD2_2_File.var"Points:1"
CD2_2_z_w = CD2_2_z/0.61

CD2_3_File = CSV.File("mesh_study/c19 cd2_3.csv")
CD2_3_h = CD2_3_File.h * 1e6
CD2_3_z = CD2_3_File.var"Points:1"
CD2_3_z_w = CD2_3_z/0.61

## Plotting

plot(CD1_1_z_w, CD1_1_h, label = "7872 cells", legend=:bottom)
plot!(CD1_2_z_w, CD1_2_h, label = "17328 cells")
plot!(CD1_3_z_w, CD1_3_h, label = "38600 cells")
xlabel!("Z/w")
ylabel!("thickness (μm)")

plot(CD2_1_z_w, CD2_1_h, label = "15050 cells", legend=:bottom)
plot!(CD2_2_z_w, CD2_2_h, label = "34528 cells")
plot!(CD2_3_z_w, CD2_3_h, label = "78450 cells")
xlabel!("Z/w")
ylabel!("thickness (μm)")

## Mesh refinement
r21_1 = (38600/17328)^(1/3)
ϕ1_1 = 0.0017901
ϕ2_1 = 0.0020089
ϕ3_1 = 0.002443
p_1 = 1/log(r21_1)*log(abs((ϕ3_1 - ϕ1_1)/(ϕ2_1 - ϕ1_1)))
ea_21_1 = abs((ϕ1_1-ϕ2_1)/ϕ1_1)
GCI_21_1 = (1.25*ea_21_1)/((r21_1)^p_1-1)

r21_2 = (78450/34528)^(1/3)
ϕ1_2 = 0.000559
ϕ2_2 = 0.0005438
ϕ3_2 = 0.0005312
p_2 = 1/log(r21_2)*log(abs((ϕ3_2 - ϕ1_2)/(ϕ2_2 - ϕ1_2)))
ea_21_2 = abs((ϕ1_2-ϕ2_2)/ϕ1_2)
GCI_21_2 = (1.25*ea_21_2)/((r21_2)^p_2-1)