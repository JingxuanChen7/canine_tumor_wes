# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:53:41 2020

The script take the depth of coverage as input file 
to identify the frequency of the coverage of depth for all of the base
position for each sample 

@author: abc73_000
"""
## look at all the positions to the end of the bases

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import skew 
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt
from matplotlib.backends.backend_pdf import PdfPages

"""
with open("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Osteo/Normal/Normal_list")as f:
    file_list = f.read().split('\n')[:-1]
    
pp = PdfPages('/Users/kun-linho/Desktop/Normal_Osteosarcoma.pdf')
"""


input_file = '/Users/kun-linho/Desktop/SRR7780741_DepthofCoverage_Distribution.txt'




#'/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Osteo/Normal/'+i
#'G:\\Pan_cancer\\Pan_cancer_mapping_result\\Distribution\\Mammary\\Normal\\SRR7780741_DepthofCoverage_Distribution.txt' 
# "/Volumes/Research_Data/Pan_cancer/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt"
#sys.argv[1]
#"/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt"
file_name = '/Users/kun-linho/Desktop/test_output.txt'

#i.split('_')[0]
#sys.argv[2]
#"SRR7780741"
#sys.argv[2]
Cancer_type = 'Mammary_Cancer'
#'Mammary_Cancer'
#sys.argv[3]
Status =  'Normal'
#'Normal'
#sys.argv[4]
## fill up the gap data, if there is no data at that position, then fill up with 0
with open(input_file,'r')as f:
    file = f.read().split('\n')[:-1]

file.remove('1 Total_Depth')
summary = {}
original_list=[]
for i in range(len(file)):
    freq = int(file[i].split(" ")[0])
    pos = int(file[i].split(" ")[1])
    summary[pos]= int(freq)
    for j in range(freq):
        original_list.append(pos)
        
total_line = len(original_list)

if original_list[-1] > 1000:
    last_pos = original_list[-1]
else:
    last_pos = 1000
    


for i in range(last_pos+1):
    if i not in summary.keys():
        summary[i]=0
        
order_summary={}
for i in sorted(summary.keys()):
    order_summary[i]=summary[i]

total_data = pd.DataFrame(list(order_summary.items()),columns=['Position', 'Frequency'])

freq_arr = total_data['Frequency'].values
pos_arr = total_data['Position'].values
#freq_arr_1000 = freq_arr[0:1000]
#pos_arr_1000 = pos_arr[0:1000]

total_array = np.array(original_list)

"""
original_list_1000=[]
for i in range(1000):
    pos = pos_arr[i]
    freq = freq_arr_1000[i]
    for j in range(freq):
        original_list_1000.append(pos)
"""

average = np.mean(np.array(original_list))
#average_1000= np.mean(np.array(original_list_1000))
std = np.std(np.array(original_list))
#std_1000 = np.std(np.array(original_list_1000))  

mu = average
#mu_1000= average_1000

prob_arr= freq_arr/total_line

#freq_list= list(freq_arr)
poisson_fract_list=[]

for i in range(last_pos+1):
    value = poisson.pmf(i,mu)
    poisson_fract_list.append(value)
    
poisson_fract_list_1000=[]

"""
for i in range(1000):
    value = poisson.pmf(i,mu_1000)
    poisson_fract_list_1000.append(value)
#[poisson.pmf(i,mu) for i in range(last_pos+1)]
"""

#pos_1000_line = len(original_list_1000)

#prob_arr_1000 = freq_arr_1000/ pos_1000_line

#Prob_data_1000 = pd.DataFrame(data= {'Prob' :prob_arr_1000,'Pos': pos_arr_1000})
pd.DataFrame(list(order_summary.items()),columns=['Position', 'Frequency'])
rmse = sqrt(mean_squared_error(prob_arr, np.array(poisson_fract_list))) ## compare the proportion
sumOfSqerror= sqrt(sum((prob_arr-np.array(poisson_fract_list))**2)) ## rmse that didn't divide the N, vs proportion

#rmse_1000 = sqrt(mean_squared_error(prob_arr_1000, np.array(poisson_fract_list_1000))) ## compare the proportion
#sumOfSqerror_1000= sqrt(sum((prob_arr_1000-np.array(poisson_fract_list_1000))**2))


poisson_list_count =[]
poisson_list_count_1000 =[]
for i in range(last_pos+1):
    value = poisson.pmf(i,mu)*total_line
    poisson_list_count.append(value)

"""
for i in range(1000):
    value = poisson.pmf(i,mu)*pos_1000_line
    poisson_list_count_1000.append(value)
"""     

rmse_count = sqrt(mean_squared_error(freq_arr, np.array(poisson_list_count)))     
sumOfSqerror_count= sqrt(sum((freq_arr-np.array(poisson_list_count))**2))

#rmse_count_1000 = sqrt(mean_squared_error(freq_arr_1000, np.array(poisson_list_count_1000)))     
#sumOfSqerror_count_1000= sqrt(sum((freq_arr_1000-np.array(poisson_list_count_1000))**2))

  
output = open(file_name+'_randomness_summary.txt','w')

output.write(file_name+'\t'+str(average)+'\t' \
+str(std)+'\t'+str(rmse)+'\t'+str(sumOfSqerror)+'\t' \
+str(Cancer_type)+'\t'+str(Status)+'\n')

output.close()


"""
p = sns.lineplot(x ='Pos', y = 'Prob', data =Prob_data_1000)
plt.title(Cancer_type+"_"+file_name+"_"+Status)
p.set(xlabel='Coverage_Depth', ylabel= "fraction of the number of position")
#sns.distplot(total_data['Frequency'],hist=False, kde=True, axlabel= 'Frequency', color='orange')
p = sns.distplot(np.array(original_list),kde=False, axlabel= 'Position', color='black', bins=200)  
plt.title(file_name)
p.set_yscale('log')   # set into log scale 
#pp.savefig()
#p = p.map(plt.hist, "value", color="r", log=True)
plt.figure(figsize=(13,6))
sns.set(font_scale=2)      
#sns.kdeplot(np.array(original_list), shade=True,bw='scott');
#pp.savefig()
plt.close() 
#pp.close()
"""    
"""   
    
    plt.figure(figsize=(13,6))
    sns.set(font_scale=2)
    #p = sns.lineplot(x ='Position', y = 'Frequency', data =total_data)
    #sns.distplot(total_data['Frequency'],hist=False, kde=True, axlabel= 'Frequency', color='orange')
    p = sns.distplot(np.array(original_list),kde=False, axlabel= 'Position', color='black', bins=200)  
    plt.title(file_name)
    p.set_yscale('log')   # set into log scale 
    pp.savefig()
    #p = p.map(plt.hist, "value", color="r", log=True)
    plt.figure(figsize=(13,6))
    sns.set(font_scale=2)      
    sns.kdeplot(np.array(original_list), shade=True,bw='scott');
    pp.savefig()
    plt.close() 
pp.close()


plt.figure(figsize=(13,6))
sns.set(font_scale=2)
p = sns.distplot(np.array(original_list_1000), hist=True,kde=False, \
             color = 'darkblue', \
             bins = 500, \
             hist_kws={'edgecolor':'black'}, \
             kde_kws={'linewidth': 3})
plt.title(file_name)
#p.set_yscale('log')   # set into log scale 
#pp.savefig()

plt.figure(figsize=(13,6))
sns.set(font_scale=2)      
sns.kdeplot(np.array(original_list_1000), shade=True,bw='scott');
## median= position (sum(freq)+1)/2 th position 
## Calcaulte Standard deviation in a same for loop ###
## sqrt(sum of (freq*(value - mean)**2)/sum of (freq)) 
"""
"""
poisson_original_list=[]
for i in range(last_pos+1):
  value =int(poisson.pmf(i,mu)*total_line)
  for j in range(value):
      poisson_original_list.append(i)



def makeOneSpace(String):
    content = String.split(' ')
    i = 0
    start = True
    while start:
        if content[1].isdigit() and content[2].isdigit():
            start = False
        else:
            i+=1
            content = String.split(' ')[i:]
            #print(i)   
    return content 
"""

