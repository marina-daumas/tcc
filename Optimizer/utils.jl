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
        # This segment generates uniform demands, of length demand_length,in the interval [-d_max,d_max] 
        if demand_type == "uniform"
            d_mat[:,j] = Int.(round.(rand(Uniform(0,dM),demand_length)))  
        end
        # This segment generates approximately normally distributed demands
        if demand_type == "normal"
            d_mat[:,j] = Int.(round.(dM .+ rand(Normal(0, std_dev), demand_length)))
        end
        # This segment generates demands with a Poisson distribution
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


function plot_data(evolution, label, ylabel, title)
    p = plot(layout = (3,2), size=(900,600), legend=:bottomright)

    nsubs = 6
    for i in 1:nsubs
        # first series for this subplot
       
        y1 = evolution[i]
        x1 = 1:length(y1)
        plot!(p[i], x1, y1, seriestype=:scatter, colour="red", label = label[i], xlabel = L"k", ylabel = ylabel[i], title = title[i])
        plot!(p[i], x1, y1, colour="red", label=false)

        # set subplot background if desired
        plot!(p[i], bg_inside = :lightblue)
    end

    display(p)
end

function plot_results(label, ylabel, title, res, metrics, i, blr, d, a)
    horiz = metrics[!, "horiz"][1]
    println(res.id)
    evolution = [res.X[i, 1:horiz,1], res.Y[i, 1:horiz,1], res.Ser[i, 1:horiz,1], blr[i, 1:horiz,1], d[1:horiz, i, 1], a[1:horiz, i, 1]]
    
    plot_data(evolution, label, ylabel, title)
end

function plot_paired_data(evolution, label, ylabel, title)
    p = plot(layout = (3,2), size=(900,600), legend=:bottomright)

    nsubs = 6
    for j in 1:nsubs
        # first series for this subplot
        i1 = j
        if i1 <= length(evolution)
            y1 = evolution[i1]
            x1 = 1:length(y1)
            plot!(p[j], x1, y1, seriestype=:scatter, colour="red", label = label[i1], xlabel = L"k", ylabel = ylabel[j], title = title[j])
            plot!(p[j], x1, y1, colour="red", label=false)
        end

        # paired series (if present) plotted on same subplot
        i2 = j + nsubs
        if i2 <= length(evolution)
            y2 = evolution[i2]
            x2 = 1:length(y2)
            plot!(p[j], x2, y2, seriestype=:scatter, colour="blue", label = label[i2])
            plot!(p[j], x2, y2, colour="blue", label=false)
        end
        
        # set subplot background if desired
        plot!(p[j], bg_inside = :lightblue)
    end

    display(p)
end

function plot_paired_results(label, ylabel, title, res, metrics, i, blr, d, a, res2, blr2, d2, a2)
    horiz = metrics[!, "horiz"][1]
    println(res.id)
    evolution = [res.X[i, 1:horiz,1], res.Y[i, 1:horiz,1], res.Ser[i, 1:horiz,1], blr[i, 1:horiz,1], d[1:horiz, i, 1], a[1:horiz, i, 1], 
        res2.X[i, 1:horiz,1], res2.Y[i, 1:horiz,1], res2.Ser[i, 1:horiz,1], blr2[i, 1:horiz,1]]
    
    plot_paired_data(evolution, label, ylabel, title)
end

function add_line_to_csv(result_id, res, horiz, d_prop, a_prop, cost, av_blr, av_dr, av_ser, total_clients, count)
    df = DataFrame(result_id=result_id, model_id=res.id, horiz=horiz, XM=res.bds.XM, YM=res.bds.YM, phiM=res.bds.phiM, serM=res.bds.serM, tserM=res.bds.tserM, 
    d_type=d_prop.type, dM=d_prop.M, d_std_dev=d_prop.std_dev, a_type=a_prop.type, aM=a_prop.M, a_std_dev=a_prop.std_dev,
    c_blr=res.c.blr, c_ser=res.c.ser, c_Z=res.c.Z, c_L=res.c.L,
    cost=cost, av_blr=av_blr, av_dr=av_dr, av_ser=av_ser, total_clients=total_clients, count=count)

    CSV.write("results.csv", df, append=true)
end


function compute_metrics(res, horiz, N)
    blr = res.L./(res.L.+res.Z)

    av_blr = sum(blr)/(N*horiz) 
    av_ser = sum(res.Ser)/(N*horiz)
    av_dr = sum(res.dr)/(N*horiz)

    total_clients = sum(res.Z[:,horiz+1])/N
    cost = res.c.ser*av_ser+ res.c.blr*av_blr

    return cost, av_blr, av_dr, av_ser, total_clients, blr
end

#### Structs ####

struct bounds
        XM
        YM
        phiM
        serM  # different for each model
        tserM # only defined for the model with delay
end

struct c
    Z
    L
    ser
    blr
end

struct initial_conditions
        X0
        Y0
        L0
        Z0
    end

struct distribution_properties
    M
    std_dev
    type
end

struct result
    id
    status_opt
    X
    Y
    Z
    L
    n
    Q
    dr
    phi
    Cin
    Ser
    Sl
    Sa
    Sst
    Sc
    Sin
    J
    c
    bds
end