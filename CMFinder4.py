from biopandas.pdb import PandasPdb
from math import *
import numpy as np
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import tqdm


def euclidean_distance(x,y):
    return sqrt(sum(pow(a-b,2) for a,b in zip(x,y)))


def getProteinData(protein):

    protein_data=PandasPdb().read_pdb(os.path.abspath(protein))

    return protein_data


def extractCoordinates(protein_data):

    protein_data=protein_data.df['ATOM']
    protein_data=protein_data[(protein_data['chain_id']=='A')]
    protein_data=protein_data[(protein_data['atom_name']=='CB')&(((protein_data['alt_loc']=='')|(protein_data['alt_loc']=='A')))|
    (protein_data['residue_name']=='GLY')&(protein_data['atom_name']=='CA')&((protein_data['alt_loc']=='')|(protein_data['alt_loc']=='A'))]
    coordinate_data=protein_data.loc[:,['x_coord','y_coord','z_coord']]
    
    return coordinate_data


def getContactmap(coordinate_data):
    
    contactmap=np.empty((100,100))
    for i in range(100):
        for j in range(100):
                residue_pair_1=coordinate_data.iloc[i][['x_coord','y_coord','z_coord']].values
                residue_pair_2=coordinate_data.iloc[j][['x_coord','y_coord','z_coord']].values
                contactmap[i][j]=euclidean_distance(residue_pair_1,residue_pair_2)
                if contactmap[i][j]<=8:
                    contactmap[i][j]=1
                else:
                    contactmap[i][j]=0

    return contactmap




def getSimilarity(contactmap1, contactmap2):
    
    similarity_list=[]

    if len(contactmap1) > len(contactmap2):
        pass
    else:
        for i in range(len(contactmap2)-len(contactmap1)+1):
            similarity_element=0
            target_element=0
            for j in range(len(contactmap1)):
                for k in range(len(contactmap1)):
                    if(contactmap1[j][k]==1):
                        target_element+=1
                        if (contactmap1[j][k]==contactmap2[j+i][k+i]):
                            similarity_element+=1
            similarity_list.append(similarity_element/target_element) # 0  

        similarity=max(similarity_list)

    return similarity
    
def main():
    proteins=sys.argv[1:]
    protein_list=[]


    for protein in proteins:
        protein_data = getProteinData(protein)
        coordinate_data = extractCoordinates(protein_data)
        contactmap = getContactmap(coordinate_data)
        protein_list.append(contactmap)


    similarity = getSimilarity(protein_list[0], protein_list[1])

    print('similarity = {:.3f}'.format(similarity))


if __name__ == '__main__':
    main()

