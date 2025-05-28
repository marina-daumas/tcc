using Statistics
using DelimitedFiles

function CC_compare_stats_nl_model(ic, bds, exd, c_ser, c_blr, demand_stats, method)
#usage: To train, set N_train \neq 0, N_test = 0, to test, vice versa.  

    # Demand characteristics(in case some demand needs to be replaced due to infeasibility)
    demand_length = demand_stats.demand_length
    dM = demand_stats.dM
    aM = demand_stats.aM
    demand_type = demand_stats.demand_type
    std_dev_a = demand_stats.std_dev_a
    std_dev_d = demand_stats.std_dev_d
    limit_count = demand_stats.limit_count # max number of trials to find feasible demands
    
    Ncols = size(exd.d_mat,2);; # number of columns of the demands matrix   

    # Optimal solutions
    X = zeros(demand_length+1, Ncols)
    Y = zeros(demand_length+1, Ncols)
    s = zeros(demand_length, Ncols)
    blr = zeros(demand_length+1, Ncols)

    # Demands loop
    for i in 1:Ncols
        println("Demand #", i)
        d = exd.d_mat[:,i]  #i'th demand from exog data matrix
        a = exd.a_mat[:,i]  #i'th abandonment vector from exog data matrix
        
        if method == "os" #one step ahead                       
            
            optimal_os = false;
            count = 1           
            
            while !optimal_os && count<=limit_count
                (optimal_os, X[:,i], Y[:,i], s[:,i], blr[:,i]) = CC_nl_fobj_os(ic,bds,c_ser,c_blr,d,a)                
            
                # If infeasibilities are found, replace demand with a new one and try again  
                if !optimal_os
                    println("At least one solution is infeasible, replacing demand ", i, "(count = ", count, "/", limit_count, ")")
                    d = exd.d_mat[:,i] = demand_stats.demand_generator_mat(1,demand_length,dM,demand_type,std_dev_d) #demand_generator_mat(3,20,14,"uniform",1)
                    a = exd.a_mat[:,i] = demand_stats.demand_generator_mat(1,demand_length,aM,demand_type,std_dev_a)
                    count = count+1;

                    if count>=limit_count
                        println("WARNING: No feasible solutions for demand ", i, " found after ", limit_count, " atempts of generating a new one")
                    end
                end              
            end

        else #omniscient
            
            optimal_om = false;
            count = 1            

            while !optimal_om && count <= limit_count
                (optimal_om, X[:,i], Y[:,i], s[:,i], blr[:,i])  = CC_nl_fobj_om(ic,bds,c_ser,c_blr,d,a);
                
                # If infeasibilities are found, replace demand with a new one and try again  
                if !optimal_om
                    println("At least one solution is infeasible, replacing demand ", i, "(count = ", count, "/", limit_count, ")")
                    d = exd.d_mat[:,i] = demand_stats.demand_generator_mat(1,demand_length,dM,demand_type,std_dev_d) 
                    a = exd.a_mat[:,i] = demand_stats.demand_generator_mat(1,demand_length,aM,demand_type,std_dev_a)                    
                    count = count+1

                    if count>=limit_count
                        println("WARNING: No feasible solutions for demand ", i, " found after ", limit_count, " atempts of generating a new one")
                    end
                end   
            end
        end        
    end

    return X, Y, s, blr
end