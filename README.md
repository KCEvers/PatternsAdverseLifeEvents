This repository contains the analysis scripts for the paper: ``Non-random patterns in the co-occurrence and accumulation of adverse life events in two national panel datasets" by Kyra Evers, Denny Borsboom, Eiko Fried, Fred Hasselman, František Bartoš, and Lourens Waldorp.

The code is organized into the following folders:
* code/analysis.R (main) runs all analyses on the SHP and HILDA data, using the following scripts:
  * code/setup.R loads the necessary libraries, sets up the parameters, and defines some helper functions
  * code/prepare_data_func.R defines functions to prepare the data from SHP and HILDA for analysis
  * cooccurrence_event_types.R runs the co-occurrence of adverse life event types analysis
  * lagged_event_counts.R runs the lagged event counts analysis
  * accumulation_event_counts.R runs the accumulation of event counts analysis
* code/descriptives.R (plots, numerical) calculates and plots descriptive statistics for the SHP and HILDA data

The resulting figures can be found in the `figs` folder.

The SHP data is available at https://www.swissubase.ch/en/catalogue/studies/6097/20179/overview, with instructions for download at https://forscenter.ch/projects/swiss-household-panel/data. Access to the HILDA dataset can be requested at https://dataverse.ada.edu.au/dataverse/DSSLongitudinalStudies. 

For any questions, please contact Kyra Evers at kyra.c.evers@gmail.com.



