
include("cc_os_original.jl")
include("cc_om_original.jl")
include("cc_os_delay.jl")
include("cc_om_delay.jl")


function optimize_list(id)

    if occursin("ori", id)
        c_ser = c_sero
        c_blr = c_blro
        bds = bdso
    else
        c_ser = c_serd
        c_blr = c_blrd
        bds = bdsd
    end

    if occursin("os", id)
        J = zeros(N, horiz) 
    else
        J = zeros(N) # cost function
    end   

    # initialize variables
    status_opt = zeros(N)    
    X = zeros(N, horiz+1)     # current number of customers in queue x(k)S
    Y = zeros(N, horiz+1)     # current number of customers in buffer y(k)  
    Z = zeros(N, horiz+1)     # custumers served  
    L = zeros(N, horiz+1)     # number of customers lost
    n = zeros(N, horiz)       # number of empty slots in the queue
    Q = zeros(N, horiz)       # custumers entering queue
    dr = zeros(N, horiz)      # dropped due to full buffer

    phi = zeros(N, horiz)       # number of customers admitted to queue
    Cin = zeros(N, horiz)       # number of customers entering server

    S = zeros(N, horiz)                   # number of active servers
    Sl = zeros(N, horiz)                  # number of free available servers
    Sa = zeros(Bool, N, horiz, bds.serM)           # server activation status (0 - active, 1 - inactive)
    Sst = zeros(Bool, N, horiz, bds.serM)          # server status (0 - free, 1 - busy)
    Sc = zeros(Bool, N, horiz+1, bds.serM, bds.tserM)  # server conveyor
    Sin = zeros(Bool, N, horiz, bds.serM, bds.tserM)   # server input
    Saux = zeros(Bool, N, horiz, bds.serM)

    for i in 1:N
        d = d_mat[1:horiz, i];  # demand for incoming calls
        a = a_mat[1:horiz, i];  # abandonment for calls

        if id == "os_ori"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], S[i,:], J[i,:] = cc_os_original(ic, bds, c_ser, c_blr, d, a, lin_fobj);
        elseif id == "om_ori"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], S[i,:], J[i] = cc_om_original(ic, bds, c_ser, c_blr, d, a, lin_fobj);
        elseif id == "os_del"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Cin[i,:], S[i,:], Sl[i,:], Sa[i,:,:], Sst[i,:,:], Sc[i,:,:,:], Sin[i,:,:,:], Saux[i,:,:], J[i,:] = cc_os_delay(ic, bds, a, d, df_input, lin_fobj);
        elseif id == "om_del"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Cin[i,:], S[i,:], Sl[i,:], Sa[i,:,:], Sst[i,:,:], Sc[i,:,:,:], Sin[i,:,:,:], Saux[i,:,:], J[i] = cc_om_delay(ic, bds, a, d, df_input, lin_fobj);
        end  
    end
    
    res = result(id, status_opt, X, Y, Z, L, n, Q, dr, phi, Cin, S, Sl, Sa, Sst, Sc, Sin, Saux, J, c_ser, c_blr);
    print(status_opt)
    return res
end