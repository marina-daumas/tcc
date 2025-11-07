include("optimizers//cc_no_delay_os.jl")
include("optimizers//cc_no_delay_om.jl")
include("optimizers//cc_delay_os.jl")
include("optimizers//cc_delay_om.jl")


include("optimizers//cc_delay_os_nonlin.jl")
include("optimizers//cc_delay_om_nonlin.jl")


function optimize_list(id, d_mat, a_mat)

    if occursin("nd", id)
        c = c_nd
        bds = bds_nd
    else
        c = c_d
        bds = bds_d
    end

    if occursin("os", id)
        J = zeros(N, horiz)
    else
        J = zeros(N) # cost function
    end

    # initialize variables
    status_opt = Bool.(zeros(N))
    X = zeros(N, horiz+1)     # current number of customers in queue x(k)S
    Y = zeros(N, horiz+1)     # current number of customers in buffer y(k)
    Z = zeros(N, horiz+1)     # custumers served
    L = zeros(N, horiz+1)     # number of customers lost
    n = zeros(N, horiz)       # number of empty slots in the queue
    Q = zeros(N, horiz)       # custumers entering queue
    dr = zeros(N, horiz)      # dropped due to full buffer

    phi = zeros(N, horiz)       # number of customers admitted to queue
    Cin = zeros(N, horiz)       # number of customers entering server

    Ser = zeros(N, horiz)                   # number of active servers
    Sl = zeros(N, horiz, bds.serM)                  # number of free available servers
    Sa = zeros(N, horiz, bds.serM)           # server activation status (0 - active, 1 - inactive)
    Sst = zeros(N, horiz, bds.serM)          # server status (0 - free, 1 - busy)
    Sc = zeros(N, horiz+1, bds.serM, bds.tserM)  # server conveyor
    Sin = zeros(N, horiz, bds.serM, bds.tserM)   # server input

    av_count = 0
    for i in 1:N

        status_opt[i] = false;
        count = 0          
            
        while !status_opt[i] && count<limit_count
            d = d_mat[1:horiz, i];  # demand for incoming calls
            a = a_mat[1:horiz, i];  # abandonment for calls
            
            if id == "nd_os" 
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Ser[i,:],J[i,:] = cc_no_delay_os(ic,bds,c,a,d);
            elseif id == "nd_om"
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Ser[i,:],J[i] = cc_no_delay_om(ic,bds,c,a,d);
            elseif id == "d_os" 
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Cin[i,:],Ser[i,:],Sl[i,:,:],Sa[i,:,:],Sst[i,:,:],Sc[i,:,:,:],Sin[i,:,:,:],J[i,:] = cc_delay_os(ic,bds,c,a,d,df_input);
            elseif id == "d_om" 
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Cin[i,:],Ser[i,:],Sl[i,:,:],Sa[i,:,:],Sst[i,:,:],Sc[i,:,:,:],Sin[i,:,:,:],J[i] = cc_delay_om(ic,bds,c,a,d,df_input);
            elseif id == "d_os_nonlin"
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Cin[i,:],Ser[i,:],Sl[i,:,:],Sa[i,:,:],Sst[i,:,:],Sc[i,:,:,:],Sin[i,:,:,:],J[i,:] = cc_delay_os_nonlin(ic,bds,c,a,d,df_input);
            elseif id == "d_om_nonlin"
                status_opt[i],X[i,:],Y[i,:],Z[i,:],L[i,:],n[i,:],Q[i,:],dr[i,:],phi[i,:],Cin[i,:],Ser[i,:],Sl[i,:,:],Sa[i,:,:],Sst[i,:,:],Sc[i,:,:,:],Sin[i,:,:,:],J[i] = cc_delay_om_nonlin(ic,bds,c,a,d,df_input);
            end                
            
            # If infeasibilities are found, replace demand with a new one and try again  
            if !status_opt[i]
                count = count+1;
                if count>=limit_count
                    println("WARNING: No feasible solutions for demand ", i, " found after ", limit_count, " atempts of generating a new one")
                else
                    d_mat[1:horiz, i] = demand_generator_mat(1, horiz, d_prop.M, d_prop.type, d_prop.std_dev); #demand_generator_mat(3,20,14,"uniform",1)
                    a_mat[1:horiz, i] = demand_generator_mat(1, horiz, a_prop.M, a_prop.type, a_prop.std_dev);
                    println("At least one solution is infeasible, replacing demand ", i, "(count = ", count, "/", limit_count, ")")
                end
                
            end              
        end
        av_count += count        
    end
    av_count = av_count/N

    res = result(id, status_opt, X, Y, Z, L, n, Q, dr, phi, Cin, Ser, Sl, Sa, Sst, Sc, Sin, J, c, bds);
    print(status_opt)
    return res, d_mat, a_mat, av_count
end