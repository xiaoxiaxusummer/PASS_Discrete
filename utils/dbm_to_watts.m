function P_watt = dbm_to_watts(P_dbm)
    P_watt = 10.^((P_dbm - 30) / 10);
end
