#Matrix demand generator: d_mat = demand_generator_mat(N_dem,demand_length,dM,demand_type,std_dev)
using Distributions
#N_dem number of demands
#dM = 15 #Maximum demand
#demand_type = "uniform", "normal" (in this case specify std_dev) or "poisson"
#This function will generate a demand matrix d_mat (dimensions matrix demand_length X N-dem), each column of which is a demand
function demand_generator_mat(N_dem,demand_length,dM,demand_type,std_dev)
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