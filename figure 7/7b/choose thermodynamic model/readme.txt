run optimizerwithlhssampling.m to generate thermodynamic models 

"workspace_gatecombinations_16ics_4.mat" generated at the end of the run contains parameter combinations that represent thermodynamic models(ground truth systems). 

since thermodynamic model for networks with high in degrees take linger to compute. user can stop the optimization as soon as four steadystates are reached denoted by objective function values being less than -40(iteratively displayed).
corresponding optimized parameter can be found in "SaveBest.mat " where all the best result for each generation is saved.

run "data_training.m" to get training data. "trainingdata.mat" will be saved at the end of the run.

run "data_for prediction.m" to generate data to be used prediction test and in plotting figure 8. 