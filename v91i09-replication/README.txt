The code for replication is divided into separate folders.

1) The analysis of gene expression levels for liver samples of female mice
   is given in the folder LiverAnalysis.

2) The analysis of HIV immuno-phenotypes is contained in the folder
   ResponderVsControllerAnalysis.

3) SimulationStudy contains the code for and the aggregated results of the
   simulations described in the article.

   The simulations were carried out using UCLA's Hoffman2 Cluster and recomputed
   at the HPC cluster at WU Vienna University of Economic and Business using
   R version 3.5.2 (detailed sessionInfo is provided in each directory). Both
   of these clusters use the SGE scheduler. The file doit.qsub specifies how the
   jobs were submitted to the scheduler via "qsub doit.qsub" on the command
   line.

   The scripts running the linear model simulations are located in the
   corresponding subdirectories. These scripts run simulations for
   fuzzy forests, conditional inference forests, and random forests.
   The results of each simulation are stored in .RData files
   of the form "out/outx" where 'x' is the index of the simulation.
   After all simulations are completed, they can be recombined via recombine.R.

   The script running the nonlinear simulation is located in 
   SimulationResults/NonLinearModelSims/nl_simspaper/, and
   SimulationResults/VimSimulations carries out simulations
   for the calculation of VIMs.
