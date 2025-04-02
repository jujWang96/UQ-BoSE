# Auto-BUQ

"Astro_sim" folder includes Auto-BUQ analysis script for NGC2300,NGC2300 Newton-XMM and Arp299 Newton-XMM data. To produce the simulation results, run "run_2300chandra_UQ_analysis.m", "run_2300XMM_UQ_analysis.m" and "run_299XMM_UQ_analysis.m". To produce plot analysis and reproduce figures, run "post_2300_UQ_analysis.m", "post_2300XMM_UQ_analysis.m" and "post_299XMM_UQ_analysis.m".

"Astro_sim/sim" folder includes analysis on merge methods on a circular disc source with four embedded point sources. Use "run_sim4.m", "run_sim4_random.m" and "run_sim4_beam.m" to generate simulation results with greedy merge, random merge and beam merge. Use "merge_result_plot_tile.m" for visualization. 

"Bulleye_sim" folder includes analysis on ESRG method. Use "run_sim_bulleye.m" to generate simulation results and "make_bulleye_plot_tile.m" for visualization.

"coverage_sim" folder includes analysis on Auto-BUQ method. Use "main_sim_coverage_tune_ellipse.m", "main_sim_coverage_tune_bridge.m" and "main_sim_coverage_tune_spiral.m" to generate simulated examples for ellipse, bridge and spiral shape, and "plot_GCR_coverage_sim.m" for visualization. To reproduce the histogram of model selection results, use "run_sim_model_selection_bridge.m" , "run_sim_model_selection_ellipse.m", "run_sim_model_selection_spiral.m" for simulation and "plot_model_selection.m" for visualization. To produce the coverage probability plot in appendix, use "run_sim_coverage_tune_spiral.m", "run_sim_coverage_tune_ellipse.m", "run_sim_coverage_tune_bridge.m" to reproduce simulated data and "plot_coverage_prob.m" for visualization. 



