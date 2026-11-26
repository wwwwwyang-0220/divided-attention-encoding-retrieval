Title: A Memory Representational Perspective on the Asymmetric Effects of Divided Attention Encoding and Retrieval

Archived on Zenodo, including the large model files: https://doi.org/10.5281/zenodo.24545658
Code repository: https://github.com/wwwwwyang-0220/divided-attention-encoding-retrieval


DATA AND ANALYSIS CODE

  Analyses_Script.R            all analyses reported in the manuscript
  Exp1_Retrieval/              Experiment 1, DA at retrieval
  Exp2_Encoding/               Experiment 2, DA at encoding
  Cross-Experiment_Analysis/   output of the comparison between the two experiments


ANALYSIS SCRIPT

Analyses_Script.R follows the order of the Results section:

  1. Effects of DA on Associative Recognition Accuracy
  2. Effects of DA on Response Latency in the Associative Recognition Task
  3. MPT analysis
  4. Simulated Replications of the Modeling Results
  5. Concurrent Arithmetic Task Performance

Within each section the analyses are run for Experiment 1, then Experiment 2,
then the cross-experiment comparison (there is no cross-experiment part for the
simulations). The data are read in at the top of the script, so after running
that part each block can be run on its own.

Fitting the MPT models and running the simulations takes a long time. By default
(refit <- FALSE at the top of the script) the saved model fits and simulation
results are loaded instead. Set refit <- TRUE to run everything from scratch.

The script uses the R packages bruceR, TreeBUGS, dplyr, openxlsx, and logspline.
TreeBUGS needs JAGS (https://mcmc-jags.sourceforge.io/). The working directory
is set to the location of the script through rstudioapi, so run it in RStudio,
or set the working directory to this folder yourself.


FILES IN THE EXPERIMENT FOLDERS

The two folders have the same structure. File names contain "DAretrieval" in
Experiment 1 and "DAencoding" in Experiment 2, written as * below.

Data

  data.xlsx
      Trial-level data of the associative recognition task (attention condition,
      probe type, accuracy, response, RT). In Experiment 1, DA is coded
      0 = full attention, 1 = addition, 2 = division, and Probe is coded
      0 = intact, 1 = related, 2 = unrelated.

  data_secondary.xlsx
      Trial-level data of the concurrent arithmetic task (accuracy and RT).

  data_mpt.xlsx
      Response frequencies per participant, in the format needed to fit the MPT
      model.

MPT model

  *_mpt.eqn
      Model equations used to fit the data.

  *_mpt_simulation.eqn
      Model equations used to fit the simulated data sets.

  mpt_*.RData
      The fitted model, including the posterior samples. These two files (about
      600 MB each) are too large for the repository. Download them from the
      v1.0 release
      (https://github.com/wwwwwyang-0220/divided-attention-encoding-retrieval/releases/tag/v1.0)
      or from Zenodo (https://doi.org/10.5281/zenodo.24545658), and put each
      one in its experiment folder.

  mpt_*.jags, mpt_*.txt
      JAGS model file and summary of the parameter estimates, both written by
      TreeBUGS.

  mpt_*_FA_Addition.xlsx, mpt_*_FA_Division.xlsx, mpt_*_Addition_Division.xlsx
      Differences in the group-level parameters between two attention conditions
      (FA = full attention).

Simulations

  ModelSimulation/Undertrue_recovery_result_001_100.rds
      100 simulated replications, with data generated from the parameter
      estimates of the fitted model.

  ModelSimulation/Undernull_recovery_result_001_100.rds
      100 simulated replications, with data generated under no difference
      between the two conditions.


CROSS-EXPERIMENT FOLDER

  Encoding_vs_Retrieval_interaction.xlsx
      For each MPT parameter, the DA effect at encoding, the DA effect at
      retrieval, and the difference between the two.
