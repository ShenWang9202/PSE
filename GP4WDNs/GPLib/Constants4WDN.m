classdef Constants4WDN
    properties( Constant = true )
        pi = 3.141592654;
        GPMperCFS= 448.831;
        AFDperCFS= 1.9837;
        MGDperCFS= 0.64632;
        IMGDperCFS=0.5382;
        LPSperCFS= 28.317;
        M2FT = 3.28084989501;
        %LPS2GMP = 15.850372483753;
        LPMperCFS= 1699.0;
        CMHperCFS= 101.94;
        CMDperCFS= 2446.6;
        MLDperCFS= 2.4466;
        M3perFT3=  0.028317;
        LperFT3=   28.317;
        MperFT=    0.3048;
        PSIperFT=  0.4333;
        KPAperPSI= 6.895;
        KWperHP=   0.7457;
        SECperDAY= 86400;
        % Index
        GPM2CFS = 448.831;%GPMperCFS;
        LPS2GMP =  15.8502313098;%448.831/28.317;
        FT2M = 0.3048;
        m2feet = 3.28084989501;%1/FT2M;
        feet2inch = 12;
        mm2m = 1000;
        mm2inch = 0.0393701;
        L2cube_m = 1000;
                
                
        TankIndex = 8;
        PumpIndex = 9;
        SpeedIndexInXX0 = 18;
        % 
        Head_Reservior = 700;
        Reservior_index = 1;
        Hp = 5;
        ReferenceHead = 838.8;
        Delta_t = 1800;%seconds 0.5hour
        SIZE_Z_K = 7;
        SIZE_X_K = 1;
        SIZE_W_K = 7;
        SIZE_U_K = 2;
        SIZE_NODE = 8;
        SIZE_D_K = 8;
        SIZE_SPEED_K = 1;
        COUNT_TANKS = 1;
    end
end

 %https://www.mathworks.com/matlabcentral/answers/358826-best-method-to-define-constants-used-by-several-functions
