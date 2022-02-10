# README for code on severe-aGVHD data 

Data was collected from 15 days before HSCT to 100 days after HSCT. 
Correspondingly, data was stored in files with 115 rows for each patient. 
To ignore pre-HSCT data, the first 15 rows, '+offset.t' would be seen frequently in the code. 


Hyper parameters setting:

Line 19: 
For our real world data we used, we set 'population = 1' for adults and 'population = 0' for children.
Here children demo data was available to play around, please set 'population = 0'

Line 29 / 34: 
For mode 'internal validation' please set 'CV = 1' and 'HO = 0';
For mode 'external validation' please set 'CV = 0' and 'HO = 1'.

When in mode CV (set 'CV = 1' and 'HO = 0'): 
Line 30: 'CV.Training.Size = 1' to get the association between 'traning size' and daGOAT's performance in CV cohort 
Line 31: 'CV.Key.Importance = 1' to get feature imporatance in CV cohort 
Line 32: 'CV.Top.Keys = 1' to get the association between the number of top N keys and daGOAT's performance in CV cohort 
Line 35: please set 'HO.Top.Keys = 0'

When in mode HO (set 'CV = 0' and 'HO = 1'): 
Line 30: please set 'CV.Training.Size = 0' 
Line 31: please set 'CV.Key.Importance = 0' 
Line 32: please set 'CV.Top.Keys = 0' 
Line 35: 'HO.Top.Keys = 1' to get the association between the number of top N keys and daGOAT's performance in holdout cohort

Line 37: 'Save.Files = 1' to save results to local machine. 


With all hyper parameters defined, you are all set to go. Just run the whole script. 

ATTENTION: 
Run the code all over again when 
mode ('internal validation' / external validation') or 
population (adults / children) shifted.
