# REDOX-PFLOTRAN
This branch is developed to investigate subsurface redox reaction control over greenhouse gas emission from Mississippi Delta wetland.
1. The initial simulation is conducted along a salinity gradient in Barataria Basin.
2. Hourly hydrological measurements to each site are used to force the reactive-transport model (Forcing folder)
3. Soil properties data for each site are used to set up the single soil column for each site (soil folder)
4. The complex redox network can be found in MRDwet_RedoxNet.py, and new reactions can be added by revising this file.
5. The reactive-transport simulation is done by run DeltaMarshGHG.py.
6. All of the other files are files that are under development for other purposes.
