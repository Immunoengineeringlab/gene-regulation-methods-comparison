To get data to run "plot_histogram_comparisons", run files in the following order.

get groundtruth data from "choose thermodynamic model" folder. The script to be executed is specified in a separate readme file in the folder.

copy the groundtruth data to data to "train Hillcube" and "train RPM" folders and execute specified scripts(specified in readme text in each folder) to get the Hillcube models and RPM models respectively, that closely approximates the groundtruth data. 

use the data files 'workspace_gatecombinations_rpm_networkA.mat' and 'workspace_gatecombinations_hillcube_networkA.mat from "train Hillcube" and "train RPM" folders to plot histogram using "plot_histogram_comparisons.m"