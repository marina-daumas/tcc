using JuMP
using HiGHS


function optimize_model(horiz, a, d, c_cin, c_ser, ic, bds, std_dev, rnd_tser)
    # bounds
    XM = bds.XM
    YM = bds.YM     
    phiM = bds.phiM
    serM = bds.serM

    J = zeros(horiz)            # cost function
    X = zeros(Int, horiz+1)     # current number of customers in queue x(k)S
    Y = zeros(Int, horiz+1)     # current number of customers in buffer y(k)  
    Z = zeros(Int, horiz+1)     # custumers served  
    L = zeros(Int, horiz+1)     # number of customers lost

    n = zeros(Int, horiz+1)     # number of empty slots in the queue
    Q = zeros(Int, horiz+1)     # custumers entering queue
    dr = zeros(Int, horiz+1)    # dropped due to full buffer

    phi = zeros(Int, horiz+1)     # number of customers admitted to queue
    Cin = zeros(Int, horiz+1)     # number of customers entering server
    Cout = zeros(Int, horiz+1)    # number of customers leaving server

    S = zeros(Int, horiz+1)                 # number of active servers
    Sa = zeros(Bool, horiz+1, serM)         # server activation status (0 - active, 1 - inactive)
    Sl = zeros(Int, horiz+1)                # number of free available servers
    Sst = zeros(Bool, horiz+1, serM)        # server status (0 - free, 1 - busy)
    Sc = zeros(Bool, horiz+1, serM, tserM)  # server conveyor
    Sin = zeros(Bool, horiz+1, serM, tserM) # server input

    Saux = zeros(Bool, horiz+1, serM)       # auxiliary server variable
    b1_opt = zeros(horiz)
    b2_opt = zeros(horiz)

    B = zeros(Float64, horiz+1)             # balking rate

    transition_matrix = zeros(Bool, tserM, tserM)
    for i in 1:tserM-1
        transition_matrix[i, i+1] = 1
    end

    if !rnd_tser
        default_input = zeros(Bool, serM, tserM)  # default input (case with constant service time)
        for i in 1:serM
            default_input[i, 1] = 1
        end
    end

    # initial conditions
    X[1] = ic.X0
    Y[1] = ic.Y0
    L[1] = ic.L0
    Z[1] = ic.Z0

    for t in 1:horiz
        cc_lin_conv = Model(HiGHS.Optimizer)
        set_silent(cc_lin_conv)

        a[t] = min(X[t], max(0, Int.(round.(X[t]/4 .+ rand(Normal(0, std_dev))))))

        default_input = zeros(Bool, serM, tserM)  
        for i in 1:serM
            slot = rand(1:tserM)  # randomly select a slot for each server
            default_input[i, slot] = 1
        end

        M1 = max(d[t], YM) + 10; # must be larger than d[t] and n[t] 
        @variable(cc_lin_conv, b1, Bin);

        M2 = max(X[t], serM) + 10; # must be larger than x[k] and Sl[k]
        @variable(cc_lin_conv, b2, Bin);

        @variable(cc_lin_conv, 0 <= XL[1:2] <= XM, Int)    
        @variable(cc_lin_conv, 0 <= YL[1:2] <= YM, Int)     
        @variable(cc_lin_conv, 0 <= ZL[1:2], Int)  
        @variable(cc_lin_conv, 0 <= LL[1:2], Int)         
                            
        @variable(cc_lin_conv, 0 <= nL <= YM, Int) 
        @variable(cc_lin_conv, 0 <= QL, Int)      
        @variable(cc_lin_conv, 0 <= drL, Int)   
        
        @variable(cc_lin_conv, 0 <= phiL <= phiM, Int)   
        @variable(cc_lin_conv, 0 <= CinL <= serM, Int)   
        @variable(cc_lin_conv, 0 <= CoutL[1:2] <= serM, Int)    

        @variable(cc_lin_conv, 0 <= SL <= serM, Int)  
        @variable(cc_lin_conv, 0 <= SaL[1:serM], Bin)
        @variable(cc_lin_conv, 0 <= SlL <= serM, Int)  
        @variable(cc_lin_conv, 0 <= SstL[1:serM], Bin)  
        @variable(cc_lin_conv, 0 <= ScL[1:2, 1:serM, 1:tserM], Bin) 
        @variable(cc_lin_conv, 0 <= SinL[1:serM, 1:tserM], Bin)
        @variable(cc_lin_conv, 0 <= SauxL[1:serM], Bin)

        # Initial conditions of buffer and queue for each optimization
        @constraint(cc_lin_conv, XL[1] == X[t])
        @constraint(cc_lin_conv, YL[1] == Y[t])
        @constraint(cc_lin_conv, ZL[1] == Z[t])
        @constraint(cc_lin_conv, LL[1] == L[t])
        @constraint(cc_lin_conv, CoutL[1] == Cout[t])
        @constraint(cc_lin_conv, ScL[1, :, :] == Sc[t, :, :])  # server conveyor status

        # Problem constraints  
        @constraint(cc_lin_conv, XL[2] == XL[1] + phiL - a[t] - CinL)
        @constraint(cc_lin_conv, YL[2] == YL[1] + QL - phiL)
        @constraint(cc_lin_conv, ZL[2] == ZL[1] + CoutL[1])
        @constraint(cc_lin_conv, LL[2] == LL[1] + drL + a[t])

        @constraint(cc_lin_conv, nL == YM - YL[1] + phiL)  # number of empty slots in the queue
        @constraint(cc_lin_conv, QL <= d[t])
        @constraint(cc_lin_conv, QL <= nL)
        @constraint(cc_lin_conv, QL >= d[t]-M1*b1)
        @constraint(cc_lin_conv, QL >= nL-(1-b1)*M1)  

        @constraint(cc_lin_conv, drL >= d[t]-nL)         
        @constraint(cc_lin_conv, drL <= d[t]-QL) 

        @constraint(cc_lin_conv, SL == serM - sum(SaL))              
        @constraint(cc_lin_conv, SstL == sum(ScL[1,:,:], dims=2))  
        @constraint(cc_lin_conv, SlL == SL - sum(SstL))    
        @constraint(cc_lin_conv, [i=1:serM, j=1:tserM], SinL[i, j] == default_input[i, j] .* SauxL[i])
        @constraint(cc_lin_conv, ScL[2, :, :] == ScL[1, :, :]*transition_matrix + SinL)

        @constraint(cc_lin_conv, CinL <= XL[1]-a[t])
        @constraint(cc_lin_conv, CinL <= SlL)
        @constraint(cc_lin_conv, CinL >= XL[1]-a[t]-M2*b2)
        @constraint(cc_lin_conv, CinL >= SlL-(1-b2)*M2)  
        @constraint(cc_lin_conv, CoutL[2] == sum(ScL[1,:,tserM]))  

        @constraint(cc_lin_conv, SauxL <= ones(Bool, serM) - SstL - SaL)  
        @constraint(cc_lin_conv, CinL == sum(SauxL))  

        # Objective function
        @objective(cc_lin_conv, Min, LL[2] + c_ser*SL - c_cin*CinL - phiL) # design better objective function

        JuMP.optimize!(cc_lin_conv)
        
        status = termination_status(cc_lin_conv)
        if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status==MOI.ALMOST_LOCALLY_SOLVED) && has_values(cc_lin_conv)
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
            Sst[t, :] = JuMP.value.(SstL);
            Sc[t+1, :, :] = JuMP.value.(ScL[2, :, :]);
            Sin[t, :, :] = JuMP.value.(SinL);
            Saux[t, :] = JuMP.value.(SauxL);

            b1_opt[t] = JuMP.value(b1);
            b2_opt[t] = JuMP.value(b2);

            J[t] = objective_value(cc_lin_conv) # Repeated line

            # Display results
            # println("Optimal solution:")
            # println("Objective value = ", objective_value(cc_lin_conv))
        else
            # println(t, " ", status)
        end
    end

    return X, Y, Z, L, n, Q, dr, phi, Cin, Cout, S, Sl, Sst, Sc, Sin, Saux, b1_opt, b2_opt, J

end