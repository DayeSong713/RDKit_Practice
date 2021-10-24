from biopandas.pdb import PandasPdb
from math import *
import numpy as np
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt


def euclidean_distance(x,y):
    return sqrt(sum(pow(a-b,2) for a,b in zip(x,y)))

if __name__ == '__main__':
    proteins=sys.argv[1:3]
    protein_list=[]


    for i in proteins:

        protein_list.append(PandasPdb().read_pdb(os.path.abspath(i)))

    

    for i in range(len(protein_list)):
       
        protein_list[i]=protein_list[i].df['ATOM']
        protein_list[i]=protein_list[i][(protein_list[i]['chain_id']=='A')]
        protein_list[i]=protein_list[i][(protein_list[i]['atom_name']=='CB')&(((protein_list[i]['alt_loc']=='')|(protein_list[i]['alt_loc']=='A')))|
        (protein_list[i]['residue_name']=='GLY')&(protein_list[i]['atom_name']=='CA')&((protein_list[i]['alt_loc']=='')|(protein_list[i]['alt_loc']=='A'))]
        protein_list[i]=protein_list[i].loc[:,['x_coord','y_coord','z_coord']]

    arr1=[]
    

    for k in range(len(protein_list)):
        arr2=np.empty((len(protein_list[k]),len(protein_list[k])))
        for i in range(len(protein_list[k])):
            for j in range(len(protein_list[k])):
                residue_pair_1=protein_list[k].iloc[i][['x_coord','y_coord','z_coord']].values
                residue_pair_2=protein_list[k].iloc[j][['x_coord','y_coord','z_coord']].values
                arr2[i][j]=euclidean_distance(residue_pair_1,residue_pair_2)
        arr1.append(arr2)


    #print(arr1)

    
    arr3=[]
    
    for k in range(len(arr1)):
        arr4=np.empty((len(arr1[k]),len(arr1[k])))
        for i in range(len(arr1[k])):
            for j in range(len(arr1[k])):
                if arr1[k][i][j]<=8:
                    arr4[i][j]=1
                else:
                    arr4[i][j]=0
        arr3.append(arr4)
    
    #print(arr3)


    similarity_list=[]



    
    for i in range(len(arr3[1][0])-len(arr3[0][0])+1):
        similarity=0
        target=0
        for j in range(len(arr3[0][0])):
            for k in range(len(arr3[0][0])):
                if(arr3[0][j][k]==1):
                    target+=1
                    if (arr3[0][j][k]==arr3[1][j+i][k+i]):
                        similarity+=1
        similarity_list.append(similarity/target) 

    print("Similarity : ",max(similarity_list))

