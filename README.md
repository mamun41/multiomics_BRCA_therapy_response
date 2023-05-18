## multiomics_BRCA_therapy_response
##Thanks all the guidance from Professor Kumar Selvarajoo and Dr Md Mamunur Rashid in Bioinformatics Institute, Astar

##Thanks all the data and code support from:
#ammut, S. J., Crispin-Ortuzar, M., Chin, S. F., Provenzano, E., Bardwell, H. A., Ma, W., ... & Caldas, C. (2022). Multi-omic machine learning predictor of breast cancer therapy response. Nature, 601(7894), 623-629.
#Singh, A., Shannon, C. P., Gautier, B., Rohart, F., Vacher, M., Tebbutt, S. J., & Lê Cao, K. A. (2019). DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. Bioinformatics, 35(17), 3055-3062.

##DIABLO: for feature selection
#4 omics data: DEG.csv, somactic mutations.csv, molecular_markers.csv, microenvironment.csv
#Patients information: Supplementary Tables.xlsx

##Machine_learning: train on features selected by DIABLO
#All training data: all.csv
#Training code: train_main.py
#Training dic: my_list***.pkl & featnames.p
#External data: external_test.csv
#Test code: run_models.py

##You should replace the file path with your own path
