
include("cc_no_delay_os.jl")
include("cc_no_delay_om.jl")
include("cc_delay_os.jl")
include("cc_delay_om.jl")

include("cc_delay_os_copy.jl")

function optimize_list(id)

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

    Ser = zeros(N, horiz)                   # number of active servers
    Sl = zeros(N, horiz)                  # number of free available servers
    Sa = zeros(Bool, N, horiz, bds.serM)           # server activation status (0 - active, 1 - inactive)
    Sst = zeros(Bool, N, horiz, bds.serM)          # server status (0 - free, 1 - busy)
    Sc = zeros(Bool, N, horiz+1, bds.serM, bds.tserM)  # server conveyor
    Sin = zeros(Bool, N, horiz, bds.serM, bds.tserM)   # server input
    Saux = zeros(Bool, N, horiz, bds.serM)

    for i in 1:N
        d = d_mat[1:horiz, i];  # demand for incoming calls
        a = a_mat[1:horiz, i];  # abandonment for calls

        if id == "nd_os"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Ser[i,:], J[i,:] = cc_no_delay_os(ic, bds, c, a, d);
        elseif id == "nd_om"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Ser[i,:], J[i] = cc_no_delay_om(ic, bds, c, a, d);
        elseif id == "d_os"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Cin[i,:], Ser[i,:], Sl[i,:], Sa[i,:,:], Sst[i,:,:], Sc[i,:,:,:], Sin[i,:,:,:], Saux[i,:,:], J[i,:] = cc_delay_os(ic, bds, c, a, d, df_input);
        elseif id == "d_om"
            status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Cin[i,:], Ser[i,:], Sl[i,:], Sa[i,:,:], Sst[i,:,:], Sc[i,:,:,:], Sin[i,:,:,:], Saux[i,:,:], J[i] = cc_delay_om(ic, bds, c, a, d, df_input);
        # else
        #     Sl = zeros(N, horiz, bds.serM)
        #     # status_opt[i], X[i,:], Y[i,:], Z[i,:], L[i,:], n[i,:], Q[i,:], dr[i,:], phi[i,:], Cin[i,:], Ser[i,:], Sl[i,:,:], Sa[i,:,:], Sst[i,:,:], Sc[i,:,:,:], Sin[i,:,:,:], J[i] = cc_delay_os_copy(ic, bds, c, a, d, df_input);
        #     status_opt[i], X[i,:], Y[i,:], Z[i,:],a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17 = cc_delay_os_copy(ic, bds, c, a, d, df_input);
        end  
    end
    
    res = result(id, status_opt, X, Y, Z, L, n, Q, dr, phi, Cin, Ser, Sl, Sa, Sst, Sc, Sin, Saux, J, c);
    print(status_opt)
    return res
end

struct result
    id
    status_opt  
    X
    Y
    Z
    L
    n
    Q
    dr
    phi
    Cin
    Ser
    Sl
    Sa
    Sst
    Sc
    Sin
    Saux
    J
    c
end