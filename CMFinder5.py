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


def getProteinData(protein_path):
    
    protein_data=PandasPdb().read_pdb(protein_path)

    return protein_data
    

def extractCoordinates(protein_data):

    protein_data=protein_data.df['ATOM']
    protein_data=protein_data[(protein_data['chain_id']=='A')]
    protein_data=protein_data[(protein_data['atom_name']=='CB')&(((protein_data['alt_loc']=='')|(protein_data['alt_loc']=='A')))|
    (protein_data['residue_name']=='GLY')&(protein_data['atom_name']=='CA')&((protein_data['alt_loc']=='')|(protein_data['alt_loc']=='A'))]
    coordinate_data=protein_data.loc[:,['x_coord','y_coord','z_coord']]

    return coordinate_data


def getContactmap(coordinate_data):
    
    contactmap=np.empty((len(coordinate_data),len(coordinate_data)))
    for i in range(len(coordinate_data)):
        for j in range(len(coordinate_data)):
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
        
        for i in range(len(contactmap1)-len(contactmap2)+1):
            similarity_element=0
            target_element=0
            for j in range(len(contactmap2)):
                for k in range(len(contactmap2)):
                    if(contactmap2[j][k]==1):
                        target_element+=1
                        if (contactmap2[j][k]==contactmap1[j+i][k+i]):
                            similarity_element+=1
            if target_element==0:
                similarity=0
            else:
                similarity=similarity_element/target_element

            similarity_list.append(similarity) # 0  

        similarity=max(similarity_list)
    
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
            if target_element==0:
                similarity=0
            else:
                similarity=similarity_element/target_element

            similarity_list.append(similarity) 

        similarity=max(similarity_list)
        
    return similarity


    
    
def main():

    
    target_protein=sys.argv[1]
    target_dir = "/home/ailon/Desktop/python_practice/PDBfile/target/"
    target_protein_path=target_dir + target_protein

    if not os.path.exists(target_protein_path): 
        print("File not found")
        exit()

    
    target_protein_data = getProteinData(target_protein_path)
    target_coordinate_data = extractCoordinates(target_protein_data)
    target_contactmap = getContactmap(target_coordinate_data)
    
    
    


    DB_dir = "/home/ailon/Desktop/python_practice/PDBfile/DB/"
    
    filelist = os.listdir(DB_dir)
    filelist = [f for f in filelist if f.endswith('.pdb')]

    # templist = []
    # for f in filelist:
    #     if f.endswith('.pdb'):
    #         templist.append(f)
    # filelist = templist
    print("similarity")
    for filename in filelist:

        db_protein_path = DB_dir + filename        
        db_protein_data = getProteinData(db_protein_path)
        db_coordinate_data = extractCoordinates(db_protein_data)
        db_contactmap = getContactmap(db_coordinate_data)
        similarity = getSimilarity(target_contactmap, db_contactmap)
        print(filename.replace(".pdb",""),": {:.3f}".format(similarity))
        

    


if __name__ == '__main__':
    main()
