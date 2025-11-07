# Call Center Model One Step Ahead (OSA) optimization with delay 
# Linear objective function

using JuMP
using HiGHS
using Ipopt

function cc_delay_os_nonlin(ic, bds, c, a, d, df_input)
    # Time horizon
    horiz = length(d)

    # Check for optimality
    optimal=false
    
    # bounds
    XM = bds.XM
    YM = bds.YM     
    phiM = bds.phiM
    serM = bds.serM
    tserM = bds.tserM

    # Output variables
    J = zeros(horiz)            # cost function
    X = zeros(horiz+1)     # current number of customers in queue x(k)S
    Y = zeros(horiz+1)     # current number of customers in buffer y(k)  
    Z = zeros(horiz+1)     # custumers served  
    L = zeros(horiz+1)     # number of customers lost

    n = zeros(horiz)     # number of empty slots in the queue
    Q = zeros(horiz)     # custumers entering queue
    dr = zeros(horiz)    # dropped due to full buffer

    phi = zeros(horiz)       # number of customers admitted to queue
    Cin = zeros(horiz)       # number of customers entering server

    S = zeros(horiz)                   # number of active servers
    Sa = zeros(horiz, serM)           # server activation status (0 - active, 1 - inactive)
    Sl = zeros(horiz, serM)                  # number of free available servers
    Sst = zeros(horiz, serM)          # server status (0 - free, 1 - busy)
    Sc = zeros(horiz+1, serM, tserM)  # server conveyor
    Sin = zeros(horiz, serM, tserM)   # server input

    blr = zeros(horiz+1) # Non lin

    transition_matrix = zeros(Bool, tserM, tserM)
    for i in 1:tserM-1
        transition_matrix[i, i+1] = 1
    end

    # initial conditions
    X[1] = ic.X0
    Y[1] = ic.Y0
    L[1] = ic.L0
    Z[1] = ic.Z0

    blr[1] = L[1]/(L[1] + Z[1]) # Non lin
    
    # Optimization loop
    for t in 1:horiz
        # Setting up solver
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 1)
        highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

        cc_os_delay = Model(
            optimizer_with_attributes(
                Juniper.Optimizer, 
                "nl_solver" => ipopt,
                "mip_solver" => highs, 
                "allow_almost_solved" => false,
                "feasibility_pump" => true        
                ))      # Non lin

        # cc_os_delay = Model(HiGHS.Optimizer)

        # No screen output
        set_silent(cc_os_delay)

        M = max(d[t], YM) + 10; # must be larger than d[t] and n[t] 
        @variable(cc_os_delay, b, Bin);

        # Optimization variables
        @variable(cc_os_delay, 0 <= XL[1:2] <= XM, Int)    
        @variable(cc_os_delay, 0 <= YL[1:2] <= YM, Int)     
        @variable(cc_os_delay, 0 <= ZL[1:2], Int)  
        @variable(cc_os_delay, 0 <= LL[1:2], Int)         
                            
        @variable(cc_os_delay, 0 <= nL <= YM, Int) 
        @variable(cc_os_delay, 0 <= QL, Int)      
        @variable(cc_os_delay, 0 <= drL, Int)   
        
        @variable(cc_os_delay, 0 <= phiL <= phiM, Int)   
        @variable(cc_os_delay, 0 <= CinL <= serM, Int)     

        @variable(cc_os_delay, 0 <= SL <= serM, Int)  
        @variable(cc_os_delay, SaL[1:serM], Bin)
        @variable(cc_os_delay, SlL[1:serM], Bin)  
        @variable(cc_os_delay, SstL[1:serM], Bin)  
        @variable(cc_os_delay, ScL[1:2, 1:serM, 1:tserM], Bin) 
        @variable(cc_os_delay, SinL[1:serM, 1:tserM], Bin)

        # Initial conditions of buffer and queue for each optimization
        @constraint(cc_os_delay, XL[1] == X[t])
        @constraint(cc_os_delay, YL[1] == Y[t])
        @constraint(cc_os_delay, ZL[1] == Z[t])
        @constraint(cc_os_delay, LL[1] == L[t])
        @constraint(cc_os_delay, ScL[1, :, :] == Sc[t, :, :])  # server conveyor status

        # Problem constraints  
        @constraint(cc_os_delay, XL[2] <= XL[1] + phiL - a[t] - CinL)
        @constraint(cc_os_delay, XL[2] >= XL[1] + phiL - a[t] - CinL)
        @constraint(cc_os_delay, YL[2] <= YL[1] + QL - phiL)
        @constraint(cc_os_delay, YL[2] >= YL[1] + QL - phiL)
        @constraint(cc_os_delay, ZL[2] <= ZL[1] + CinL)
        @constraint(cc_os_delay, ZL[2] >= ZL[1] + CinL)
        @constraint(cc_os_delay, LL[2] <= LL[1] + drL + a[t])
        @constraint(cc_os_delay, LL[2] >= LL[1] + drL + a[t])
        

        @constraint(cc_os_delay, nL == YM - YL[1] + phiL)  # number of empty slots in the queue
        @constraint(cc_os_delay, QL <= d[t])
        @constraint(cc_os_delay, QL <= nL)
        @constraint(cc_os_delay, QL >= d[t]-M*b)
        @constraint(cc_os_delay, QL >= nL-(1-b)*M)  

        @constraint(cc_os_delay, drL >= d[t]-nL)         
        @constraint(cc_os_delay, drL <= d[t]-QL) 

        @constraint(cc_os_delay, SL <= serM - sum(SaL))     
        @constraint(cc_os_delay, SL >= serM - sum(SaL))  

        @constraint(cc_os_delay, SstL .<= sum(ScL[1,:,:], dims=2))  
        @constraint(cc_os_delay, SstL .>= sum(ScL[1,:,:], dims=2))
        @constraint(cc_os_delay, SlL <= ones(Bool, serM) - SstL - SaL)
        @constraint(cc_os_delay, ScL[2, :, :] .>= ScL[1, :, :]*transition_matrix + SinL)

        @constraint(cc_os_delay, ScL[2, :, :] .<= ScL[1, :, :]*transition_matrix + SinL)
        @constraint(cc_os_delay, [i=1:serM, j=1:tserM], SinL[i, j] <= df_input[t, i, j] .* SlL[i])
        @constraint(cc_os_delay, [i=1:serM, j=1:tserM], SinL[i, j] >= df_input[t, i, j] .* SlL[i])

        @constraint(cc_os_delay, CinL <= sum(SlL))        
        
        # Objective function
        # lin_fobj(S, dr, Cin, phi, Z, L) = c.ser*sum(S) - c.Z*sum(Z) +c.L*sum(L)
        # @expression(cc_os_delay, expr, lin_fobj(SL, drL, CinL, phiL, ZL[2], LL[2])) 
        # @objective(cc_os_delay, Min, expr)

        @NLobjective(cc_os_delay, Min, c.ser*SL + c.blr*LL[2]/(LL[2] + ZL[2])) # Non lin
        
        # Compute solution  
        JuMP.optimize!(cc_os_delay)
        
        # Check if solver found optimal solution
        status = termination_status(cc_os_delay)
        if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_os_delay)
            # Save computed values
            X[t+1] = round.(JuMP.value(XL[2]));
            Y[t+1] = round.(JuMP.value(YL[2]));
            Z[t+1] = round.(JuMP.value(ZL[2]));
            L[t+1] = round.(JuMP.value(LL[2]));

            n[t] = round.(JuMP.value(nL));
            Q[t] = round.(JuMP.value(QL));
            dr[t] = round.(JuMP.value(drL));

            phi[t] = round.(JuMP.value(phiL));
            Cin[t] = round.(JuMP.value(CinL));

            S[t] = round.(JuMP.value(SL));
            Sl[t,:] = round.(JuMP.value.(SlL));
            Sa[t,:] = round.(JuMP.value.(SaL));
            Sst[t,:] = round.(JuMP.value.(SstL));
            Sc[t+1,:,:] = round.(JuMP.value.(ScL[2, :, :]));
            Sin[t,:,:] = round.(JuMP.value.(SinL));

            J[t] = objective_value(cc_os_delay)

            blr[t+1] = L[t+1]/(L[t+1] + Z[t+1]) # Non lin
            optimal = true;
        else
             # Informs the calling function the problem is infeasible for this demand
            optimal = false;
            println(status)
            break
        end
    end
    return optimal, X, Y, Z, L, n, Q, dr, phi, Cin, S, Sl, Sa, Sst, Sc, Sin, J

end