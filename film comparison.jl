using CSV
using Plots

file_my_data = CSV.File("poor simulation data.csv")
my_h_Data = file_my_data.h
h_m = file_my_data.h.*1e6
my_z_Pos = file_my_data.var"Points:1"
z_w = (my_z_Pos.+0.305)./0.610

file_paper_data = CSV.File("plot-data.csv")
paper_h_Data = file_paper_data.var" y"
paper_z_w = file_paper_data.x

plot(z_w, h_m, label="my data")
plot!(paper_z_w, paper_h_Data, label="paper data")
xlabel!("z/W")
ylabel!("Film Height (μm)")