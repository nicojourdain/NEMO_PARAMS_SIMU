## Jourdain et al. (2022)

Paper submitted to GRL

* NEMO-3.6 : dev\_r5151\_UKMO\_ISF branch at revision 5932 (see fortran files saved in ../Planchat\_2022/WORK)

* Very few personal customizations (see ../Planchat\_2022/MY\_SRC).

* Configuration name: AMUXL12.L75 (1/12°, 75 vertical levels).

* Produced by Nicolas Jourdain (IGE, Grenoble, France) on occigen (CINES) in 2020 and 2021 (HPC allocation by GENCI).

* Original simulation names: 
   - nemo\_AMUXL12\_GNJ002\_BM02MAR (1979-2009, ensemble member A)
   - nemo\_AMUXL12\_GNJ002\_BM03MAR (1979-2009, ensemble member B)
   - nemo\_AMUXL12\_GNJ002\_BM04MAR (1979-2009, ensemble member C)
   - nemo\_AMUXL12\_GNJ002\_BM02MARrcp85 ("2070-2100", ensemble member A)
   - nemo\_AMUXL12\_GNJ002\_BM03MARrcBDY ("2070-2100", ensemble member B)
   - nemo\_AMUXL12\_GNJ002\_BM03MARrcp85 ("2070-2100", ensemble member B, no BDY perturbation)
   - nemo\_AMUXL12\_GNJ002\_BM04MARrcp85 ("2070-2100", ensemble member C)
   - nemo\_AMUXL12\_GNJ002\_BM04MARrcICB ("2070-2100", ensemble member C, increased iceberg melt)

* See ocean and sea-ice model parameters in files:
   - cpp\_AMU12r\_nodiaharm.fcm
   - namelist\_ice\_nemo\_GENERIC\_AMUXL12\_A (ensemble member A)
   - namelist\_ice\_nemo\_GENERIC\_AMUXL12\_B (ensemble member B)
   - namelist\_ice\_nemo\_GENERIC\_AMUXL12\_C (ensemble member C)
   - namelist\_nemo\_GENERIC\_AMUXL12\_A (ensemble member A)
   - namelist\_nemo\_GENERIC\_AMUXL12\_B (ensemble member B)
   - namelist\_nemo\_GENERIC\_AMUXL12\_C (ensemble member C)

* Output and forcing files available on request to <nicolas.jourdain@univ-grenoble-alpes.fr> 
