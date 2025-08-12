function S2R_Settings()

    par_in = Dict(

    "special_implementation" => "swissrad", # "none" or "swissrad" 
        
    # transmissivit/radiation settings
    "calc_trans" => true,
    "calc_swr"   => 0, # 0 = off; 1 = potential swr (atm_trans = 1)

    "t1"    => "01.06.2020 00:00:00", # "dd.mm.yyyy HH:MM:SS"
    "t2"    => "30.06.2020 23:40:00",
    "tstep" => 10,

    "time_zone"   => 1,
    "coor_system" => "CH1903+",

    # shis for calculation
    "SHI_summer"    => true,
    "SHI_winter"    => false,
    "SHI_terrain"   => true,
    "SHI_evergreen" => false,

    # run settings
    "batch"        => true # running in parallel or single process

    )

    # dat_in = Dict(

    # )

    return par_in


end
