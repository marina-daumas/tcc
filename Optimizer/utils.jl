using DelimitedFiles
using PrettyTables
using Distributions
using LaTeXStrings
using Plots
using CSV
using DataFrames


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


function df_input_generator(horiz, serM, tserM)
    final_input = zeros(horiz, serM, tserM)
    for i in 1:horiz
        rnd_input = zeros(serM, tserM)
        rnd_index = Int.(round.(rand(Uniform(1,tserM), serM)))
        for j in 1:serM
            slot = rnd_index[j]  
            rnd_input[j, slot] = 1
        end
        final_input[i,:,:] = rnd_input
    end

    return final_input
end


function printTable(data, header)
    println("Results")
    pretty_table(data; header=header,formatters=ft_printf("%5.3f",1:11))
end


function plotData(evolution, ylabel, title)
    gr()
    p = plot(evolution, layout = 4, seriestype=:scatter, label=false,xlabel=[L"k" L"k" L"k" L"k"],ylabel=ylabel, title=title,
        palette=cgrad.([:grays :blues :heat :lightrainbow]), bg_inside=[:lightblue :lightblue :lightblue :lightblue])
    display(p)
end


function plot_results(res, horiz, i)
    println(res.id)
    evolution = [res.X[i, 1:horiz,1], res.Z[i, 1:horiz,1], res.Ser[i, 1:horiz,1], res.L[i, 1:horiz,1]]
    ylabel = [L"X" L"Z" L"S" L"L"]
    title = ["queue occupancy" "total customers served " "active servers" "lost customers"]
    plotData(evolution, ylabel, title)
end


function add_line_to_csv(result_id, res, horiz, d_prop, a_prop, cost, av_blr, av_dr, av_ser, total_clients)
    df = DataFrame(result_id=result_id, model_id=res.id, horiz=horiz, XM=res.bds.XM, YM=res.bds.YM, phiM=res.bds.phiM, serM=res.bds.serM, tserM=res.bds.tserM, 
    d_type=d_prop.type, dM=d_prop.M, d_std_dev=d_prop.std_dev, a_type=a_prop.type, aM=a_prop.M, a_std_dev=a_prop.std_dev,
    c_blr=res.c.blr, c_ser=res.c.ser, c_Z=res.c.Z, c_L=res.c.L,
    cost=cost, av_blr=av_blr, av_dr=av_dr, av_ser=av_ser, total_clients=total_clients)

    CSV.write("results.csv", df, append=true)
end


function compute_metrics(res, horiz, N)
    blr = res.L./(res.L.+res.Z)

    av_blr = sum(blr)/(N*horiz) 
    av_ser = sum(res.Ser)/(N*horiz)
    av_dr = sum(res.dr)/(N*horiz)

    total_clients = sum(res.Z[:,horiz+1])/N
    cost = res.c.ser*av_ser+ res.c.blr*av_blr

    return cost, av_blr, av_dr, av_ser, total_clients
end
