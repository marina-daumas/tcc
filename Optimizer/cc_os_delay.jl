using JuMP
using HiGHS


function optimize_cc_os_delay_1(ic, bds, a, d, df_input, fobj)
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

    phi = zeros(Int, horiz)       # number of customers admitted to queue
    Cin = zeros(Int, horiz)       # number of customers entering server
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

    for t in 1:horiz
        cc_os_delay = Model(HiGHS.Optimizer)
        set_silent(cc_os_delay)

        M1 = max(d[t], YM) + 10; # must be larger than d[t] and n[t] 
        @variable(cc_os_delay, b1, Bin);

        # M2 = max(X[t], serM) + 10; # must be larger than x[k] and Sl[k]
        # @variable(cc_os_delay, b2, Bin);

        @variable(cc_os_delay, 0 <= XL[1:2] <= XM, Int)    
        @variable(cc_os_delay, 0 <= YL[1:2] <= YM, Int)     
        @variable(cc_os_delay, 0 <= ZL[1:2], Int)  
        @variable(cc_os_delay, 0 <= LL[1:2], Int)         
                            
        @variable(cc_os_delay, 0 <= nL <= YM, Int) 
        @variable(cc_os_delay, 0 <= QL, Int)      
        @variable(cc_os_delay, 0 <= drL, Int)   
        
        @variable(cc_os_delay, 0 <= phiL <= phiM, Int)   
        @variable(cc_os_delay, 0 <= CinL <= serM, Int)   
        @variable(cc_os_delay, 0 <= CoutL[1:2] <= serM, Int)    

        @variable(cc_os_delay, 0 <= SL <= serM, Int)  
        @variable(cc_os_delay, 0 <= SaL[1:serM], Bin)
        @variable(cc_os_delay, 0 <= SlL <= serM, Int)  
        @variable(cc_os_delay, 0 <= SstL[1:serM], Bin)  
        @variable(cc_os_delay, 0 <= ScL[1:2, 1:serM, 1:tserM], Bin) 
        @variable(cc_os_delay, 0 <= SinL[1:serM, 1:tserM], Bin)
        @variable(cc_os_delay, 0 <= SauxL[1:serM], Bin)

        # Initial conditions of buffer and queue for each optimization
        @constraint(cc_os_delay, XL[1] == X[t])
        @constraint(cc_os_delay, YL[1] == Y[t])
        @constraint(cc_os_delay, ZL[1] == Z[t])
        @constraint(cc_os_delay, LL[1] == L[t])
        @constraint(cc_os_delay, CoutL[1] == Cout[t])
        @constraint(cc_os_delay, ScL[1, :, :] == Sc[t, :, :])  # server conveyor status

        # Problem constraints  
        @constraint(cc_os_delay, XL[2] == XL[1] + phiL - a[t] - CinL)
        @constraint(cc_os_delay, YL[2] == YL[1] + QL - phiL)
        @constraint(cc_os_delay, ZL[2] == ZL[1] + CinL)
        @constraint(cc_os_delay, LL[2] == LL[1] + drL + a[t])

        @constraint(cc_os_delay, nL == YM - YL[1] + phiL)  # number of empty slots in the queue
        @constraint(cc_os_delay, QL <= d[t])
        @constraint(cc_os_delay, QL <= nL)
        @constraint(cc_os_delay, QL >= d[t]-M1*b1)
        @constraint(cc_os_delay, QL >= nL-(1-b1)*M1)  

        @constraint(cc_os_delay, drL >= d[t]-nL)         
        @constraint(cc_os_delay, drL <= d[t]-QL) 

        @constraint(cc_os_delay, SL == serM - sum(SaL))              
        @constraint(cc_os_delay, SstL == sum(ScL[1,:,:], dims=2))  
        @constraint(cc_os_delay, SlL == SL - sum(SstL))    
        @constraint(cc_os_delay, [i=1:serM, j=1:tserM], SinL[i, j] == df_input[t, i, j] .* SauxL[i])
        @constraint(cc_os_delay, ScL[2, :, :] == ScL[1, :, :]*transition_matrix + SinL)

        # @constraint(cc_os_delay, CinL <= XL[1]-a[t])
        @constraint(cc_os_delay, CinL <= SlL)
        # @constraint(cc_os_delay, CinL >= XL[1]-a[t]-M2*b2)
        # @constraint(cc_os_delay, CinL >= SlL-(1-b2)*M2)  
        @constraint(cc_os_delay, CoutL[2] == sum(ScL[1,:,tserM]))  

        @constraint(cc_os_delay, SauxL <= ones(Bool, serM) - SstL - SaL)  
        @constraint(cc_os_delay, CinL == sum(SauxL))  

        # Objective function
        @expression(cc_os_delay, expr, fobj(SL, drL, CinL, phiL, ZL[2], LL[2])) 
        @objective(cc_os_delay, Min, expr)

        JuMP.optimize!(cc_os_delay)
        
        status = termination_status(cc_os_delay)
        if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_os_delay)
            # Save computed values
            X[t+1] = JuMP.value(XL[2]);
            Y[t+1] = JuMP.value(YL[2]);
            Z[t+1] = JuMP.value(ZL[2]);
            L[t+1] = JuMP.value(LL[2]);

            n[t] = JuMP.value(nL);
            Q[t] = JuMP.value(QL);  
            dr[t] = JuMP.value(drL);
            
            phi[t] = JuMP.value(phiL);
            Cin[t] = JuMP.value(CinL);
            Cout[t+1] = JuMP.value(CoutL[2]);

            S[t] = JuMP.value(SL);
            Sl[t] = JuMP.value(SlL);
            Sa[t] = JuMP.value(SaL);
            Sst[t, :] = JuMP.value.(SstL);
            Sc[t+1, :, :] = JuMP.value.(ScL[2, :, :]);
            Sin[t, :, :] = JuMP.value.(SinL);
            Saux[t, :] = JuMP.value.(SauxL);

            b1_opt[t] = JuMP.value(b1);
            b2_opt[t] = JuMP.value(b2);

            J[t] = objective_value(cc_os_delay)

            optimal = true;
            # Display results
            # println("Optimal solution:")
            # println("Objective value = ", objective_value(cc_os_delay))
        else
            optimal = false;
            println(t, " ", status)
            break
        end
    end

    return optimal, X, Y, Z, L, n, Q, dr, phi, Cin, Cout, S, Sl, Sa, Sst, Sc, Sin, Saux, b1_opt, b2_opt, J

end