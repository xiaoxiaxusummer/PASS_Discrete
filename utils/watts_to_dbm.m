function P_dbm = watts_to_dbm(P_watt)
    P_dbm = log10(P_watt)*10+30;
end
