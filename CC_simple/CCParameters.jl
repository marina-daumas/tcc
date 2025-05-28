module CCParameters
    # Objective function coefficients
    c_blr = 140;
    c_ser = 1;
    
    # Upper bounds on X,Y,phi,ser
    XM = 5 #5
    YM =  10 #15
    phiM = 9 #15
    serM = 5 #10

    #Initial conditions
    Y0 = 2
    X0 = 4
    L0 = 0
    Z0 = 1

    struct initial_conditions
        X0
        Y0
        L0
        Z0
    end
    ic = initial_conditions(X0,Y0,L0,Z0)
    
    struct bounds
        XM
        YM
        phiM
        serM
    end
    bds = bounds(XM,YM,phiM,serM)

    struct exog_data
        d_mat
        a_mat
    end    
    exd = exog_data(Main.d_mat,Main.a_mat);    
 export ic, bds, exd, c_ser
end