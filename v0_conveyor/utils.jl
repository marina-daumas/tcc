using DelimitedFiles
using PrettyTables
using Distributions
using LaTeXStrings
using Plots

function demand_generator_mat(N_dem, demand_length, dM, demand_type, std_dev)
    d_mat = zeros(demand_length, N_dem)
    for j in 1:N_dem
        #This segment generates uniform demands, of length demand_length,in the interval [-d_max,d_max] 
        if demand_type == "uniform"
            d_mat[:,j] = Int.(round.(rand(Uniform(0,dM),demand_length)))  
        end
        ##This segment generates approximately normally distributed demands
        if demand_type == "normal"
            d_mat[:,j] = Int.(round.(dM .+ rand(Normal(0, std_dev), demand_length)))
        end
        ##This segment generates demands with a Poisson distribution
        if demand_type == "poisson"
            d_mat[:,j] = Int.(round.(dM .+ rand(Poisson(0, std_dev), demand_length)))
         end
    end
    return d_mat           
end


function printTable(data, header)
    println("Results")
    pretty_table(data; header=header,formatters=ft_printf("%5.3f",1:11))
end


function plotData(evolution, ylabel, title)
    gr()
    plot(evolution, layout = 4, seriestype=:scatter, label=false,xlabel=[L"k" L"k" L"k" L"k"],ylabel=ylabel, title=title,
        palette=cgrad.([:grays :blues :heat :lightrainbow]), bg_inside=[:lightblue :lightblue :lightblue :lightblue])
end