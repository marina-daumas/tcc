# Call Center Model One Step Ahead (OSA) optimization
# Nonlinear objective function

using JuMP
using Ipopt
using HiGHS
using Juniper

function cc_os_original(ic,bds,c_ser,c_blr,d,a,nl,fobj)
    # Time horizon
    horiz = length(d)

    # Check for optimality
    optimal=false

    # Output variables
    J = zeros(horiz)     #Cost function
    X  = zeros(horiz+1)  #Server queue length, XM is max server queue length
    Y = zeros(horiz+1)   #Buffer queue length, YM is max buffer capacity
    ser = zeros(horiz)   #Number of servers, ser
    phi = zeros(horiz)   #Transfer rate from buffer to server queue
    dr = zeros(horiz)    #Dropping rate of new arrivals to buffer
    n  = zeros(horiz)    #Number of empty spots in buffer queue
    Q  = zeros(horiz)    #Transfer rate from demand into input buffer
    L  = zeros(horiz+1)  #Lost (dropped)
    Z  = zeros(horiz+1)  #Served
    b_opt = zeros(horiz) #Binary variable for big-M logic to algebra trick
    blr = zeros(horiz+1)


    # Initial values of display variables
    X[1] = ic.X0
    Y[1] = ic.Y0
    L[1] = ic.L0
    Z[1] = ic.Z0
    blr[1] = L[1]/(L[1] + Z[1])
    
    # Optimization loop
    for t in 1:horiz  
 
        # Setting up solvers and their attributes
        # Optimization problem

        if nl
            ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 1)
            highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
            cc_nl_fobj_os = Model(
                optimizer_with_attributes(
                    Juniper.Optimizer, 
                    "nl_solver" => ipopt,
                    "mip_solver" => highs, 
                    "allow_almost_solved" => false,
                    "feasibility_pump" => true        
                    )
            )      
        else
            cc_nl_fobj_os = Model(HiGHS.Optimizer)
        end
        
        # No screen output
        set_silent(cc_nl_fobj_os)
        
        # For the integer constraint on variable b
        M = max(d[t], bds.YM)+10; # must be larger than d[t] and n[t]    
        
        # Optimization variables
        @variable(cc_nl_fobj_os, 0 <= XL[1:2] <= bds.XM, Int)
        @variable(cc_nl_fobj_os, 0 <= YL[1:2] <= bds.YM, Int)
        @variable(cc_nl_fobj_os, 0 <= serL <= bds.serM, Int)
        @variable(cc_nl_fobj_os, 0 <= phiL <= bds.phiM, Int)
        @variable(cc_nl_fobj_os, 0 <= nL <= bds.YM, Int)
        @variable(cc_nl_fobj_os, 0 <= drL, Int)
        @variable(cc_nl_fobj_os, 0 <= QL, Int)
        @variable(cc_nl_fobj_os, b , Bin);
        @variable(cc_nl_fobj_os, 0<=ZL[1:2], Int)
        @variable(cc_nl_fobj_os, 0<=LL[1:2], Int)

        # Initial conditions of buffer and queue for each optimization
        @constraint(cc_nl_fobj_os, XL[1] == X[t])
        @constraint(cc_nl_fobj_os, YL[1] == Y[t])
        @constraint(cc_nl_fobj_os, LL[1] == L[t])
        @constraint(cc_nl_fobj_os, ZL[1] == Z[t])

        # Problem constraints    
        @constraint(cc_nl_fobj_os, XL[2] == XL[1] + phiL - serL - a[t] )
        @constraint(cc_nl_fobj_os, YL[2] == YL[1] + QL - phiL )
        @constraint(cc_nl_fobj_os, nL   == bds.YM + phiL - YL[1])      
        @constraint(cc_nl_fobj_os, drL >= d[t]-nL)         
        @constraint(cc_nl_fobj_os, drL <= d[t]-QL) 

        @constraint(cc_nl_fobj_os, QL <= d[t])
        @constraint(cc_nl_fobj_os, QL <= nL)
        @constraint(cc_nl_fobj_os, QL >= d[t]-M*b )
        @constraint(cc_nl_fobj_os, QL >= nL-(1-b)*M)    

        @constraint(cc_nl_fobj_os, LL[2] == LL[1] + drL + a[t])
        @constraint(cc_nl_fobj_os, ZL[2] == ZL[1] + serL)         

        @constraint(cc_nl_fobj_os, 0 <= serL <= bds.serM)
        #@constraint(CallCenter_L_osa, serM - 2 <= serL <= serM)
        @constraint(cc_nl_fobj_os, 0 <= phiL <= bds.phiM)
        @constraint(cc_nl_fobj_os, 0 <= XL[2] <= bds.XM)
        @constraint(cc_nl_fobj_os, 0 <= YL[2] <= bds.YM)
       
        # Objective function     
        if nl   
            @NLobjective(cc_nl_fobj_os, Min, c_ser*serL + c_blr*LL[2]/(LL[2] + ZL[2]) )
        else
            @expression(cc_nl_fobj_os, expr, fobj(serL, drL, serL, phiL, ZL[2], LL[2])) 
            @objective(cc_nl_fobj_os, Min, expr)  
        end

        # Compute solution  
        JuMP.optimize!(cc_nl_fobj_os)
        
        # Check if solver found optimal solution
        status = termination_status(cc_nl_fobj_os)
        if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status==MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_nl_fobj_os)
            # Get and save solution of each optimization problem
            X[t+1] = JuMP.value(XL[2]);
            Y[t+1] = JuMP.value(YL[2]);
            L[t+1] = JuMP.value(LL[2]);
            Z[t+1] = JuMP.value(ZL[2]);
            blr[t+1] = L[t+1]/(L[t+1] + Z[t+1])
            dr[t]  = JuMP.value(drL);
            n[t]   = JuMP.value(nL);
            ser[t] = JuMP.value(serL);
            phi[t] = JuMP.value(phiL);
            Q[t]   = JuMP.value(QL);  
            b_opt[t] = JuMP.value(b);

            J[t] =  objective_value(cc_nl_fobj_os)

            # Informs the calling function an optimal solution was found
            optimal = true;               
        else
            println("Problem is infeasible...")
            
            X = zeros(horiz+1)
            Y = zeros(horiz+1)
            ser = zeros(horiz)
            blr = zeros(horiz+1)
    
            # Informs the calling function the problem is infeasible for this demand
            optimal = false;
            break
        end
    end
    print("ok")

    return optimal, X, Y, Z, L, n, Q, dr, phi, ser, b_opt, J 
end