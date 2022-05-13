using JLD2
using DataFrames, XLSX, CSV
using Plots

Ghia = DataFrame(XLSX.readtable("Re1000_ux.xlsx", "Sheet1")...)
FEM = DataFrame(XLSX.readtable("Re1000_ux_num.xlsx", "Sheet1")...)

plot(title="Lid Driven Cavity Flow, Re=1000", Ghia.y, Ghia.ux, seriestype = :scatter, label="Ghia results")
xlabel!("x[m]")
ylabel!("ux[m/s]")
plot!(FEM.arc_length,FEM.uhx, linestyle=:dot, label="FEM")
plot!(legend=:bottomright)
savefig("LidDriven_Re_1000.pdf")