using AdsorptionModel

col_params = COL_PARAMS_YOUNG();
sorb_params = SORB_PARAMS_LEWATIT();

output = simulate_process(; col_params, sorb_params, # specify column and sorbent parameters
                            duration_adsorption=8000,   duration_desorption=20000, duration_heating=2400, 
                            T_amb_adsorption  = 288.15,    T_feed_adsorption = 288.15,
                            T_amb_desorption  = 373.15,    T_feed_desorption = 373.15, 
                            u_feed_adsorption = 7.06e-2,   u_feed_desorption = 0.0, 
                            P_out_adsorption  = 1.01325e5, P_out_desorption  = 0.2e5,
                            y_H2O_adsorption  = 0.00935,   y_H2O_desorption  = 0.0,
                            steady_state_tol  = 1e-6,      max_cycles = 4, # Override default to run 4 cycles exactly
                            enable_logging=true);