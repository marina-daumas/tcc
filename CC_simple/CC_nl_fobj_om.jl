# Call Center Model Omniscient optimization
# Nonlinear objective function

using JuMP
using Ipopt
using HiGHS
using Juniper
using GLPK #not being used 

function CC_nl_fobj_om(ic,bds,c_ser,c_blr,d,a)
    # Check for optimality
    optimal = false

    # Output variables
    horiz = length(d)
    X_om  = zeros(horiz+1)
    Y_om = zeros(horiz+1)
    ser_om = zeros(horiz)
    phi_om = zeros(horiz)
    dr_om = zeros(horiz)
    n_om  = zeros(horiz)
    Q_om  = zeros(horiz)
    L_om  = zeros(horiz+1)
    Z_om  = zeros(horiz+1)
    J_om = 0 #zeros(horiz)
    blr_om = zeros(horiz+1)

    #BLK_om = zeros(horiz+1)
    # Initial values of display variables
    X_om[1] = ic.X0
    Y_om[1] = ic.Y0
    L_om[1] = ic.L0
    Z_om[1] = ic.Z0    

    #BLK_om[1] = (L0/(L0+Z0))
    # Optimization problem
    # Setting up solvers and their attributes
     # Optimization problem
     ipopt  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
     mipSolver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
     #mipSolver = optimizer_with_attributes(GLPK.Optimizer)
     
     cc_nl_fobj_om = Model(
         optimizer_with_attributes(
             Juniper.Optimizer, 
             "nl_solver" => ipopt,
             "mip_solver" => mipSolver, 
             "allow_almost_solved" => false,
             "feasibility_pump" => true
             )
     )
         
    #=
    call_center_L_om_v1 = Model(HiGHS.Optimizer)
    set_optimizer_attribute(call_center_L_om_v1, "run_crossover", "on")
    set_optimizer_attribute(call_center_L_om_v1, "presolve", "off")
    =#
    # No screen output
    set_silent(cc_nl_fobj_om)
    
    # Objective function
    fobj(s,L,Z) = c_ser*sum(s) + c_blr*sum( L./(L .+ Z) )
       
    # Optimization variables
    #=
    @variable(cc_nl_fobj_om, 0 <= XL_om[1:horiz+1] <= bds.XM, Int)
    @variable(cc_nl_fobj_om, 0 <= YL_om[1:horiz+1] <= bds.YM, Int)
    @variable(cc_nl_fobj_om, 0 <= serL_om[1:horiz] <= bds.serM, Int) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= phiL_om[1:horiz] <= bds.phiM, Int) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= nL_om[1:horiz] <= bds.YM, Int)
    @variable(cc_nl_fobj_om, 0 <= drL_om[1:horiz], Int)
    @variable(cc_nl_fobj_om, 0 <= QL_om[1:horiz], Int)
    #@variable(call_center_L_om_v1, 0<= b[1:horiz] <= 1, Int)
    @variable(cc_nl_fobj_om, b[1:horiz], Bin)
    @variable(cc_nl_fobj_om, 0<=ZL_om[1:horiz+1], Int)
    @variable(cc_nl_fobj_om, 0<=LL_om[1:horiz+1], Int)
    =#
    @variable(cc_nl_fobj_om, 0 <= XL_om[1:horiz+1] <= bds.XM, Int)
    @variable(cc_nl_fobj_om, 0 <= YL_om[1:horiz+1] <= bds.YM, Int)
    @variable(cc_nl_fobj_om, 0 <= serL_om[1:horiz] <= bds.serM) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= phiL_om[1:horiz] <= bds.phiM) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= nL_om[1:horiz] <= bds.YM)
    @variable(cc_nl_fobj_om, 0 <= drL_om[1:horiz])
    @variable(cc_nl_fobj_om, 0 <= QL_om[1:horiz])
    #@variable(call_center_L_om_v1, 0<= b[1:horiz] <= 1, Int)
    @variable(cc_nl_fobj_om, b[1:horiz], Bin)
    @variable(cc_nl_fobj_om, 0<=ZL_om[1:horiz+1])
    @variable(cc_nl_fobj_om, 0<=LL_om[1:horiz+1])    

    # Initial conditions
    @constraint(cc_nl_fobj_om, XL_om[1] == X_om[1])
    @constraint(cc_nl_fobj_om, YL_om[1] == Y_om[1])     
    @constraint(cc_nl_fobj_om, LL_om[1] == L_om[1])
    @constraint(cc_nl_fobj_om, ZL_om[1] == Z_om[1])

    # Optimization loop
    for t in 1:horiz  
        # For the integer constraint on variable b
        M = max(d[t], bds.YM)+10; # must be larger than d[t] and n[t]    
        # Problem constraints    
        @constraint(cc_nl_fobj_om, XL_om[t+1] == XL_om[t] + phiL_om[t] - serL_om[t] - a[t])
        @constraint(cc_nl_fobj_om, YL_om[t+1] == YL_om[t] + QL_om[t] - phiL_om[t] )     
        @constraint(cc_nl_fobj_om, nL_om[t]   == bds.YM + phiL_om[t] - YL_om[t])      
        @constraint(cc_nl_fobj_om, drL_om[t] >= d[t]-nL_om[t])         
        @constraint(cc_nl_fobj_om, drL_om[t] <= d[t] - QL_om[t]) 
        #@constraint(call_center_L_om_v1, drL_om[t] <= exd.d[t] )
        @constraint(cc_nl_fobj_om, QL_om[t] <= d[t])
        @constraint(cc_nl_fobj_om, QL_om[t] <= nL_om[t])
        @constraint(cc_nl_fobj_om, QL_om[t] >= d[t]-M*b[t] )
        @constraint(cc_nl_fobj_om, QL_om[t] >= nL_om[t]-(1-b[t])*M)    

        @constraint(cc_nl_fobj_om, LL_om[t+1] == LL_om[t] + drL_om[t] + a[t])
        @constraint(cc_nl_fobj_om, ZL_om[t+1] == ZL_om[t] + serL_om[t])       

        @constraint(cc_nl_fobj_om, 0 <= serL_om[t] <= bds.serM)
        @constraint(cc_nl_fobj_om, 0 <= phiL_om[t] <= bds.phiM)
        @constraint(cc_nl_fobj_om, 0 <= XL_om[t+1] <= bds.XM)
        @constraint(cc_nl_fobj_om, 0 <= YL_om[t+1] <= bds.YM)
    end
    
    # Encapsulate objective function in a @expression
    @expression(cc_nl_fobj_om, expr, fobj(serL_om, LL_om, ZL_om)  )    
    @objective(cc_nl_fobj_om, Min, expr) # uses the expression to define the nonlinear objective function    
    #@objective(cc_nl_fobj_om, Min, c_ser*sum(serL_om) + c_blr*sum(LL_om./(LL_om .+ ZL_om))) # uses the expression to define the nonlinear objective function    
    # Compute solution    
    JuMP.optimize!(cc_nl_fobj_om)
    
    # Print the model
    #print(cc_nl_fobj_om)
    
    # Save objective value
    J_om = objective_value(cc_nl_fobj_om)
 
# Check termination status
    #println("#### t = " * string(t) *" ###")

    #print("Termination status: ")
    #println(termination_status(cc_nl_fobj_om))
    
    status = termination_status(cc_nl_fobj_om)
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status==MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_nl_fobj_om)
        # Get and save solution of each optimization problem
        
        X_om = JuMP.value.(XL_om);
        Y_om = JuMP.value.(YL_om);
        ser_om = JuMP.value.(serL_om);
        L_om = JuMP.value.(LL_om);
        Z_om = JuMP.value.(ZL_om);
        blr_om = L_om ./(L_om .+ Z_om);
        optimal=true
    else
        X_om = zeros(horiz+1)
        Y_om = zeros(horiz+1)
        ser_om = zeros(horiz)
        blr_om = zeros(horiz+1)
        optimal = false
    end
    #X_om = value.(XL_om);
    #Y_om = value.(YL_om);
    #L_om = JuMP.value.(LL_om);
    #Z_om = JuMP.value.(ZL_om);
    #dr_om = value.(drL_om);
    #n_om = value.(nL_om);
    #ser_om = value.(serL_om);
    #phi_om = value.(phiL_om);
    #Q_om = value.(QL_om);
    #b_om = value.(b)
   # blr_om = L_om./(L_om.+ Z_om)
    #J_om = cpr.c_x*(X_om[1:length(ser_om)]) + cpr.c_y*(Y_om[1:length(ser_om)]) + cpr.c_ser*(ser_om) + cpr.c_dr*(dr_om)
    #omd = om_data(n_om, dr_om, b_om, Q_om, Y_om, phi_om, X_om, ser_om, blr_om, J_om)

    return optimal, X_om, Y_om, ser_om, blr_om
end