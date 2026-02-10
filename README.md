# GeneSHAP
Computing asymmetric Shapley values for high-dimensional features

This repository contains R-scripts to reproduce the results in the manuscript: "How important are the genes to explain the outcome - 
the asymmetric Shapley value as an honest importance metric for high-dimensional features" by Mark van de Wiel, Jeroen Goedhart, 
Martin Jullum, Kjersti Aas.

A. Source functions: GeneShap_source.R

B. Data example on colorectal cancer survival data.
  1. DataExample1.R (Data preprocessing + model fitting for three models)
  2. DataExample2.R (Computation of Shapley values for blockForest model)
  3. DataExample3.R (SAGE + Inference)
  4. DataExampleAll.R (Invokes all three modules)
  5. 500_CRC.Rdata (Data)

C. Simulation results, low-dimensional: SimExample_lowdim.R

D. Importance sampling: ImpSampling_check.R

