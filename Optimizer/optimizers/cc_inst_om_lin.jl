# Call Center Model Omniscient optimization
# Nonlinear objective function

using JuMP
using Ipopt
using HiGHS
using Juniper
using GLPK #not being used 

function cc_inst_om_lin(ic, bds, c, a, d)
    # Time horizon
    horiz = length(d)

    # Check for optimality
    optimal=false

    # bounds
    XM = bds.XM
    YM = bds.YM     
    phiM = bds.phiM
    serM = bds.serM

    # Output variables
    horiz = length(d)
    X  = zeros(horiz+1)
    Y = zeros(horiz+1)
    ser = zeros(horiz)
    phi = zeros(horiz)
    dr = zeros(horiz)
    n  = zeros(horiz)
    Q  = zeros(horiz)
    L  = zeros(horiz+1)
    Z  = zeros(horiz+1)
    J = 0 
    blr = zeros(horiz+1)
    b = zeros(horiz)

    # Initial values of display variables
    X[1] = ic.X0
    Y[1] = ic.Y0
    L[1] = ic.L0
    Z[1] = ic.Z0    
    
    # Setting up solvers and their attributes
    cc_nl_fobj_om = Model(HiGHS.Optimizer)
  
    # No screen output
    set_silent(cc_nl_fobj_om)
       
    @variable(cc_nl_fobj_om, 0 <= XL[1:horiz+1] <= XM, Int)
    @variable(cc_nl_fobj_om, 0 <= YL[1:horiz+1] <= YM, Int)
    @variable(cc_nl_fobj_om, 0 <= serL[1:horiz] <= serM) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= phiL[1:horiz] <= phiM) #This variable was declared as an integer
    @variable(cc_nl_fobj_om, 0 <= nL[1:horiz] <= YM)
    @variable(cc_nl_fobj_om, 0 <= drL[1:horiz])
    @variable(cc_nl_fobj_om, 0 <= QL[1:horiz])

    @variable(cc_nl_fobj_om, b[1:horiz], Bin)
    @variable(cc_nl_fobj_om, 0<=ZL[1:horiz+1])
    @variable(cc_nl_fobj_om, 0<=LL[1:horiz+1])    

    # Initial conditions
    @constraint(cc_nl_fobj_om, XL[1] == X[1])
    @constraint(cc_nl_fobj_om, YL[1] == Y[1])     
    @constraint(cc_nl_fobj_om, LL[1] == L[1])
    @constraint(cc_nl_fobj_om, ZL[1] == Z[1])

    # Optimization loop
    for t in 1:horiz  
        # For the integer constraint on variable b
        M = max(d[t], YM)+10; # must be larger than d[t] and n[t]    
        # Problem constraints    
        @constraint(cc_nl_fobj_om, XL[t+1] == XL[t] + phiL[t] - serL[t] - a[t])
        @constraint(cc_nl_fobj_om, YL[t+1] == YL[t] + QL[t] - phiL[t] )     
        @constraint(cc_nl_fobj_om, nL[t]   == YM + phiL[t] - YL[t])      
        @constraint(cc_nl_fobj_om, drL[t] >= d[t] - nL[t])         
        @constraint(cc_nl_fobj_om, drL[t] <= d[t] - QL[t]) 
        
        @constraint(cc_nl_fobj_om, QL[t] <= d[t])
        @constraint(cc_nl_fobj_om, QL[t] <= nL[t])
        @constraint(cc_nl_fobj_om, QL[t] >= d[t]-M*b[t] )
        @constraint(cc_nl_fobj_om, QL[t] >= nL[t]-(1-b[t])*M)    

        @constraint(cc_nl_fobj_om, LL[t+1] == LL[t] + drL[t] + a[t])
        @constraint(cc_nl_fobj_om, ZL[t+1] == ZL[t] + serL[t])       

        @constraint(cc_nl_fobj_om, 0 <= serL[t] <= serM)
        @constraint(cc_nl_fobj_om, 0 <= phiL[t] <= phiM)
        @constraint(cc_nl_fobj_om, 0 <= XL[t+1] <= XM)
        @constraint(cc_nl_fobj_om, 0 <= YL[t+1] <= YM)
    end
    

    lin_fobj(S, dr, phi, Z, L) = c.ser*sum(S) - c.Z*sum(Z) +c.L*sum(L)
    @expression(cc_nl_fobj_om, expr, lin_fobj(serL, drL, phiL, ZL, LL)) 
    @objective(cc_nl_fobj_om, Min, expr)
    
    # Compute solution    
    JuMP.optimize!(cc_nl_fobj_om)
    
    # Save objective value
    J = objective_value(cc_nl_fobj_om)
    
    status = termination_status(cc_nl_fobj_om)
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status==MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_nl_fobj_om)
        # Get and save solution of each optimization problem

        X = round.(JuMP.value.(XL));
        Y = round.(JuMP.value.(YL));
        ser = round.(JuMP.value.(serL));
        L = round.(JuMP.value.(LL));
        Z = round.(JuMP.value.(ZL));

        phi = round.(JuMP.value.(phiL));
        n = round.(JuMP.value.(nL));
        Q = round.(JuMP.value.(QL));
        b = round.(JuMP.value.(b));
        dr = round.(JuMP.value.(drL));

        blr = L ./(L .+ Z);

        optimal = true
    else
        optimal = false
        println(status)
    end
   
    return optimal, X, Y, Z, L, n, Q, dr, phi, ser, J
end