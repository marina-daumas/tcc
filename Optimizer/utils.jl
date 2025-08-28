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


function df_input_generator(horiz, serM, tserM, rnd_tser)
    final_input = zeros(horiz, serM, tserM)
    default_input = zeros(serM, tserM)

    for i in 1:serM
        default_input[i, 1] = 1
    end

    KK = zeros(horiz, serM)
    for i in 1:horiz
        if rnd_tser
            rnd_input = zeros(serM, tserM)
            rnd_index = Int.(round.(rand(Uniform(1,tserM), serM)))
            for j in 1:serM
                slot = rnd_index[j]  
                rnd_input[j, slot] = 1
            end
            final_input[i,:,:] = rnd_input
        else
            final_input[i,:,:] = default_input
        end
    end

    return final_input
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


function plot_results(res, horiz, i)
    println(res.id)
    evolution = [res.X[i, 1:horiz,1], res.Z[i, 1:horiz,1], res.S[i, 1:horiz,1], res.L[i, 1:horiz,1]]
    ylabel = [L"X" L"Z" L"S" L"L"]
    title = ["queue occupancy" "total customers" "active servers" "lost customers"]
    plotData(evolution, ylabel, title)
end