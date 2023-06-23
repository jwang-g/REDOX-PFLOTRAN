1. 2V3166Results, inlcude inhibition of other substrates for methanenogenesis
2. 3V3166Results, exclude inhibition of other substrates for methanenogenesis (3-yr run, 2011-2012)
3. 03-07-2022, update deltamarsh_profile to add elevated sea level senario to compare with current status, and results is saved in 2V3166Results.
4. 03-09-2022, update deltamarsh_profile to add elevated sea level senario of 0.1m
5. 03-14-2022, instead of changing with inundation depth, set fix advection rate of 3.9e-10 under wet condition, which one order higher than dry condition
6. 03-15-2022, increase salinity to 5ppt for two stations
7. 03-18-2022, set salinity to the obervation level and update temperature dependance of microbial reactions and fix a bug on rateconstants(rateconstants was overwritten somehow and was very small)
8. 03-18-2022, rerun with rateconstant bug fixed and withou bug fixed and results is named as Run_Methane__and RRun_Methane; they have the same results, so it is not a bug
9. 03-23-2022, run without oxygen inhibition term for Hydrolysis and with Songjie's water CH4 for wet boundary condition; and result is saved as Methane_output-03-23.
10. 03-23-2022, (Save as ATM_Methane_output) add oxygen inhibition term with a -1d10 threshold, because if no oxygen inhibition, DOM is high at surface;in this run, the gas species in water are equlibrium with atm, instead of using measured values.
11. 03-24-2022, (save as SH_Methane_output) using Songjie's CH4 and HCO3 water data as wet boundary condition and run with oxygen inhibition term of -1d10 threshold.
12. 03-25-2022, (save as Temp_Methane) add Q10=2.0 and using atm-water equilibrium wet boundary
13. 03-27-2022(removed), (save as wnd_Methane) add Q10=2.0 and wind effect when it is inundated and using a small CH4 background concentration from literatures which is much smaller than Songjie's data and higher than equilibrium because GoM has natural seeps
14. 03-27-2022(removed), failed; (save as wnd2_Methane) add Q10=2.0 and wind effect for exposure and inundation and using a small CH4 background concentration from literatures which is much smaller than Songjie's data and higher than equilibrium because GoM has natural seeps
15. 03-28-2022, removed; (save as wnd_Methane) add Q10=2.0 and wind effect for exposure and inundation and using a small CH4 background concentration from literatures which is much smaller than Songjie's data and higher than equilibrium because GoM has natural seeps
16. 03-29-2022, removed; (save as wnd_Methane) add Q10=2.0 and wind effect is the same for exposure and inundation and using a small CH4 background concentration from literatures which is much smaller than Songjie's data and higher than equilibrium because GoM has natural seeps
17. 03-29-2022, (save as wnd_Methane) add Q10=2.0,temperature effect on diffusion and wind effect is only for inundation and using a small CH4 background concentration from literatures which is much smaller than Songjie's data and higher than equilibrium because GoM has natural seeps
18. 03-30-2022, (save as notemp_wnd_Methane) do not use temperature effect but with wind effect and temperature effect on diffusion is commented out.
19. 03-30-2022, (save as wnd_temp_Methane) use temperature effect and with wind effect and temperature effect on diffusion is commented out.
20. 03-30-2022, (save as sal_temp_Methane) increase salinity to 5ppt for two stations, use temperature effect and wind effect and no temperature effect on diffusion.
21. 03-30-2022, (save as flood_temp_Methane) elevated sea level senario of 0.1m for two stations, use temperature effect and wind effect and no temperature effect on diffusion.
22. 03-31-2022, (save as flood_temp_Methane) elevated sea level senario of 0.1m for two stations, use temperature effect and wind effect and no temperature effect on diffusion and comment out oxy inhibition of hydrolysis.
23. 04-01-2022, (save as flood_temp_Methane) elevated sea level senario of 0.1m for two stations, use temperature effect and wind effect and no temperature effect on diffusion and with positive oxy inhibition of hydrolysis.
24. 04-01-2022, (save as sal_temp_Methane) increase salinity to 5ppt for two stations, use temperature effect and wind effect and no temperature effect on diffusion and with positive oxy inhibition of hydrolysis.
25. 04-01-2022, (save as wnd_temp_Methane) use temperature effect and wind effect and no temperature effect on diffusion and with positive oxy inhibition of hydrolysis.
26. 4-05-2022, (save as notemp_Methane) no temperature effect and wind effect and no temperature effect on diffusion and with positive oxy inhibition of hydrolysis (3166 threshold changed from 1d-11 to 1d-12 because 1d-11 didn't work).
27. 04-06-2022, (as as wnd_tmp for all sites) temperature dependency and wind effect but with negative hydrolysis oxygen inhibition, because positive or no oxygen do not work for all stations
28. 04-16-2022, (use networkv3) without inhibition terms of Fe3+, NO3- and SO4-- for methane production (as as wnd_tmp for all sites in 3V3166Results folder) temperature dependency and wind effect but with negative hydrolysis oxygen inhibition, because positive or no oxygen do not work for all stations
29. 04-17-2022, senario runs with (flood (0.1m), salinity increase(5 ppt) )(use networkv3) without inhibition terms of Fe3+, NO3- and SO4-- for methane production (as as wnd_tmp for all sites in 3V3166Results folder) temperature dependency and wind effect but with negative hydrolysis oxygen inhibition, because positive or no oxygen do not work for all stations
30. 04-18-2022, senario run with (both flood (0.1m) and salinity increase(5 ppt) )(use networkv3) without inhibition terms of Fe3+, NO3- and SO4-- for methane production (as as wnd_tmp for all sites in 3V3166Results folder) temperature dependency and wind effect but with negative hydrolysis oxygen inhibition, because positive or no oxygen do not work for all stations
## revision for JAMES journal

31. 06-05-2023, (this run is under scenarios with 0.1m sea level rise and 5ppt salinity increase) correct networkv3 reactions and keep other information the same as run 04-16-2022 (hanford.dat has been updated with more primary species)
32. 06-06-2023, (current sea level and salinity data) correct networkv3 reactions and keep other information the same as run 04-16-2022 (hanford.dat has been updated with more primary species) and run is saved in 5V3166Results.