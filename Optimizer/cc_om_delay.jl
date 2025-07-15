using JuMP
using HiGHS


function optimize_cc_om_delay(ic, bds, a, d, df_input, fobj)
    horiz = length(d)  
    optimal = false

    # bounds
    XM = bds.XM
    YM = bds.YM     
    phiM = bds.phiM
    serM = bds.serM
    tserM = bds.tserM

    J = zeros(horiz)            # cost function
    X = zeros(Int, horiz+1)     # current number of customers in queue x(k)S
    Y = zeros(Int, horiz+1)     # current number of customers in buffer y(k)  
    Z = zeros(Int, horiz+1)     # custumers served  
    L = zeros(Int, horiz+1)     # number of customers lost

    n = zeros(Int, horiz)     # number of empty slots in the queue
    Q = zeros(Int, horiz)     # custumers entering queue
    dr = zeros(Int, horiz)    # dropped due to full buffer

    phi = zeros(Int, horiz)     # number of customers admitted to queue
    Cin = zeros(Int, horiz)     # number of customers entering server
    Cout = zeros(Int, horiz+1)    # number of customers leaving server

    S = zeros(Int, horiz)                   # number of active servers
    Sa = zeros(Bool, horiz, serM)           # server activation status (0 - active, 1 - inactive)
    Sl = zeros(Int, horiz)                  # number of free available servers
    Sst = zeros(Bool, horiz, serM)          # server status (0 - free, 1 - busy)
    Sc = zeros(Bool, horiz+1, serM, tserM)  # server conveyor
    Sin = zeros(Bool, horiz, serM, tserM)   # server input

    Saux = zeros(Bool, horiz, serM)       # auxiliary server variable
    b1_opt = zeros(horiz)
    # b2_opt = zeros(horiz)

    transition_matrix = zeros(Bool, tserM, tserM)
    for i in 1:tserM-1
        transition_matrix[i, i+1] = 1
    end

    # initial conditions
    X[1] = ic.X0
    Y[1] = ic.Y0
    L[1] = ic.L0
    Z[1] = ic.Z0

    cc_om_delay = Model(HiGHS.Optimizer)
    set_silent(cc_om_delay)

    @variable(cc_om_delay, b1[1:horiz], Bin);
    @variable(cc_om_delay, b2[1:horiz], Bin);

    @variable(cc_om_delay, 0 <= XL[1:horiz+1] <= XM, Int)    
    @variable(cc_om_delay, 0 <= YL[1:horiz+1] <= YM, Int)     
    @variable(cc_om_delay, 0 <= ZL[1:horiz+1], Int)  
    @variable(cc_om_delay, 0 <= LL[1:horiz+1], Int)         
                            
    @variable(cc_om_delay, 0 <= nL[1:horiz] <= YM, Int) 
    @variable(cc_om_delay, 0 <= QL[1:horiz], Int)      
    @variable(cc_om_delay, 0 <= drL[1:horiz], Int)   
        
    @variable(cc_om_delay, 0 <= phiL[1:horiz] <= phiM, Int)   
    @variable(cc_om_delay, 0 <= CinL[1:horiz] <= serM, Int)   
    @variable(cc_om_delay, 0 <= CoutL[1:horiz+1] <= serM, Int)    

    @variable(cc_om_delay, 0 <= SL[1:horiz] <= serM, Int)  
    @variable(cc_om_delay, 0 <= SaL[1:horiz, 1:serM], Bin)
    @variable(cc_om_delay, 0 <= SlL[1:horiz] <= serM, Int)  
    @variable(cc_om_delay, 0 <= SstL[1:horiz, 1:serM], Bin)  
    @variable(cc_om_delay, 0 <= ScL[1:horiz+1, 1:serM, 1:tserM], Bin) 
    @variable(cc_om_delay, 0 <= SinL[1:horiz, 1:serM, 1:tserM], Bin)
    @variable(cc_om_delay, 0 <= SauxL[1:horiz, 1:serM], Bin)

    @constraint(cc_om_delay, XL[1] == X[1])
    @constraint(cc_om_delay, YL[1] == Y[1])
    @constraint(cc_om_delay, ZL[1] == Z[1])
    @constraint(cc_om_delay, LL[1] == L[1])

    for t in 1:horiz
        M1 = max(d[t], YM) + 10; # must be larger than d[t] and n[t] 
        # M2 = max(X[t], serM) + 10; # must be larger than x[k] and Sl[k]

        # Problem constraints  
        @constraint(cc_om_delay, XL[t+1] == XL[t] + phiL[t] - a[t] - CinL[t])
        @constraint(cc_om_delay, YL[t+1] == YL[t] + QL[t] - phiL[t])
        @constraint(cc_om_delay, ZL[t+1] == ZL[t] + CinL[t])
        @constraint(cc_om_delay, LL[t+1] == LL[t] + drL[t] + a[t])

        @constraint(cc_om_delay, nL[t] == YM - YL[t] + phiL[t])  # number of empty slots in the queue
        @constraint(cc_om_delay, QL[t] <= d[t])
        @constraint(cc_om_delay, QL[t] <= nL[t])
        @constraint(cc_om_delay, QL[t] >= d[t]-M1*b1[t])
        @constraint(cc_om_delay, QL[t] >= nL[t]-(1-b1[t])*M1)  

        @constraint(cc_om_delay, drL[t] >= d[t]-nL[t])         
        @constraint(cc_om_delay, drL[t] <= d[t]-QL[t]) 

        @constraint(cc_om_delay, SL[t] == serM - sum(SaL[t,:]))              
        @constraint(cc_om_delay, SstL[t, :] == sum(ScL[t,:,:], dims=2))  
        @constraint(cc_om_delay, SlL[t] == SL[t] - sum(SstL[t,:]))    
        @constraint(cc_om_delay, [i=1:serM, j=1:tserM], SinL[t, i, j] == df_input[t, i, j] .* SauxL[t, i])
        @constraint(cc_om_delay, ScL[t+1,:,:] == ScL[t,:,:]*transition_matrix + SinL[t,:,:])

        @constraint(cc_om_delay, CinL[t] <= XL[t]-a[t])
        @constraint(cc_om_delay, CinL[t] <= SlL[t])
        @constraint(cc_om_delay, CinL[t] >= XL[t]-a[t]-M2*b2[t])
        @constraint(cc_om_delay, CinL[t] >= SlL[t]-(1-b2[t])*M2)  
        @constraint(cc_om_delay, CoutL[t+1] == sum(ScL[t,:,tserM]))  

        @constraint(cc_om_delay, SauxL[t,:] <= ones(Bool, serM) - SstL[t,:] - SaL[t,:])  
        @constraint(cc_om_delay, CinL[t] == sum(SauxL[t,:]))  
    end

    # Objective function
    @expression(cc_om_delay, expr, fobj(SL, drL, CinL, phiL, ZL, LL)) 
    @objective(cc_om_delay, Min, expr)

    JuMP.optimize!(cc_om_delay)

    print("ok")
    status = termination_status(cc_om_delay)
    if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status==MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_om_delay)
        # Save computed values
        X = JuMP.value.(XL);
        Y = JuMP.value.(YL);
        Z = JuMP.value.(ZL);
        L = JuMP.value.(LL);

        n = JuMP.value.(nL);
        Q = JuMP.value.(QL);  
        dr = JuMP.value.(drL);
        
        phi = JuMP.value.(phiL);
        Cin = JuMP.value.(CinL);
        Cout = JuMP.value.(CoutL);

        S = JuMP.value.(SL);
        Sl = JuMP.value.(SlL);
        Sa = JuMP.value.(SaL);
        Sst = JuMP.value.(SstL);
        Sc = JuMP.value.(ScL);
        Sin = JuMP.value.(SinL);
        Saux = JuMP.value.(SauxL);

        b1_opt = JuMP.value.(b1);
        b2_opt = JuMP.value.(b2);

        J = objective_value(cc_om_delay)
        
        optimal = true;
        # Display results
        # println("Optimal solution:")
        # println("Objective value = ", objective_value(cc_om_delay))
    else        
        optimal = false
        println(status)
    end

    return optimal, X, Y, Z, L, n, Q, dr, phi, Cin, Cout, S, Sl, Sa, Sst, Sc, Sin, Saux, b1_opt, b2_opt, J

end