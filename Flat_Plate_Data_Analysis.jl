using CSV
using Plots

file = CSV.File(open("test.csv"))

velocity=10

delta_star
for i = 1:length(file.var"Points:1")
    
    if (file.var"U:0"[i] > velocity)
        global delta_star = file.var"Points:1"[i]
        break
    end
end

plot((file.var"U:0"/velocity), (file.var"Points:1"/delta_star))
plot!(xlabel="u/U∞", ylabel="y/δ")

reference = CSV.File(open("External_dataset_2.csv"), header=false)
plot!(reference.var"Column1", reference.var"Column2")