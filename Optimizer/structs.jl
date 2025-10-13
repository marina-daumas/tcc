struct bounds
        XM
        YM
        phiM
        serM  # different for each model
        tserM # only defined for the model with delay
end

struct c
    Z
    L
    ser
    blr
end

struct initial_conditions
        X0
        Y0
        L0
        Z0
    end

struct distribution_properties
    M
    std_dev
    type
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
    J
    c
    bds
end