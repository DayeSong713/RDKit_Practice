from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import matplotlib.pyplot as plt

#데이터 입력

Drug_list=['Axitinib','Lenvatinib','Alectinib','Cobimetinib','Osimertinib','Erdafitinib','Acalabrutinib','Ixazomib','Glasdegib','Palbociclib','Rucaparib','Copanlisib','Cariprazine','Rolapitant','Daclatasvir','Isavuconazonium_sulfate','Cholic_acid','Ivabradine','Deoxycholic_acid','Eluxadoline']

SMILES_list=['CNC(=O)C1=C(SC2=CC=C3C(NN=C3\C=C\C3=CC=CC=N3)=C2)C=CC=C1','COC1=C(C=C2C(OC3=CC=C(NC(=O)NC4CC4)C(Cl)=C3)=CC=NC2=C1)C(N)=O','CCC1=CC2=C(C=C1N1CCC(CC1)N1CCOCC1)C(C)(C)C1=C(C3=CC=C(C=C3N1)C#N)C2=O','OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@@H]1CCCCN1','COC1=C(NC2=NC=CC(=N2)C2=CN(C)C3=C2C=CC=C3)C=C(NC(=O)C=C)C(=C1)N(C)CCN(C)C','COC1=CC(=CC(OC)=C1)N(CCNC(C)C)C1=CC=C2N=CC(=NC2=C1)C1=CN(C)N=C1','CC#CC(=O)N1CCC[C@H]1C1=NC(=C2N1C=CN=C2N)C1=CC=C(C=C1)C(=O)NC1=CC=CC=N1','CC(C)C[C@H](NC(=O)CNC(=O)C1=CC(Cl)=CC=C1Cl)B(O)O','CN1CC[C@H](C[C@@H]1C1=NC2=CC=CC=C2N1)NC(=O)NC1=CC=C(C=C1)C#N','CC(=O)C1=C(C)C2=CN=C(NC3=NC=C(C=C3)N3CCNCC3)N=C2N(C2CCCC2)C1=O','CNCC1=CC=C(C=C1)C1=C2CCNC(=O)C3=C2C(N1)=CC(F)=C3','COC1=C(OCCCN2CCOCC2)C=CC2=C1N=C(NC(=O)C1=CN=C(N)N=C1)N1CCN=C21', 'CN(C)C(=O)N[C@H]1CC[C@H](CCN2CCN(CC2)C2=C(Cl)C(Cl)=CC=C2)CC1','C[C@@H](OC[C@]1(CC[C@]2(CCC(=O)N2)CN1)C1=CC=CC=C1)C1=CC(=CC(=C1)C(F)(F)F)C(F)(F)F','COC(=O)N[C@@H](C(C)C)C(=O)N1CCC[C@H]1C1=NC=C(N1)C1=CC=C(C=C1)C1=CC=C(C=C1)C1=CN=C(N1)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)OC)C(C)C', '[H]C(C)(OC(=O)N(C)C1=C(COC(=O)CNC)C=CC=N1)[N+]1=CN(C[C@](O)(C2=C(F)C=CC(F)=C2)[C@@]([H])(C)C2=NC(=CS2)C2=CC=C(C=C2)C#N)N=C1','[H][C@@](C)(CCC(O)=O)[C@@]1([H])CC[C@@]2([H])[C@]3([H])[C@]([H])(O)C[C@]4([H])C[C@]([H])(O)CC[C@]4(C)[C@@]3([H])C[C@]([H])(O)[C@]12C','COC1=C(OC)C=C2[C@@H](CN(C)CCCN3CCC4=CC(OC)=C(OC)C=C4CC3=O)CC2=C1','[H][C@@]12CC[C@H]([C@H](C)CCC(O)=O)[C@@]1(C)[C@@H](O)C[C@@]1([H])[C@@]2([H])CC[C@]2([H])C[C@H](O)CC[C@]12C', 'COC1=CC=C(CN([C@@H](C)C2=NC(=CN2)C2=CC=CC=C2)C(=O)[C@@H](N)CC2=C(C)C=C(C=C2C)C(N)=O)C=C1C(O)=O']

#SMILES 변환

drugs=list(map(Chem.MolFromSmiles,SMILES_list))

#Fingerprint 값 계산

Fingerprints=[FingerprintMols.FingerprintMol(x) for x in drugs]


#Tanimoto Similarity 계산 후 전체 데이터 차트 그리기

similarity_list=[]


for i in range (len(Drug_list)):
    for j in range(i+1,len(Drug_list)):
        Fingerprint_1=Fingerprints[i]
        Fingerprint_2=Fingerprints[j]
        similarity_list.append(round(DataStructs.FingerprintSimilarity(Fingerprint_1,Fingerprint_2),2))

plt.hist(similarity_list)
plt.xlabel('Tanimoto_similarity')
plt.xlim([0,1.0])
plt.show()

#Tanimoto similarity 계산 후 개별 데이터 차트 그리기

for i in range (len(Drug_list)):
    similarity_list2=[]
    for j in range (len(Drug_list)):
        Fingerprint_1=Fingerprints[i]
        Fingerprint_2=Fingerprints[j]
        similarity_list2.append(round(DataStructs.FingerprintSimilarity(Fingerprint_1,Fingerprint_2),2))
    similarity_list2.remove(1.0)
    plt.subplot(4,5,i+1)
    plt.hist(similarity_list2)
    plt.xlim([0,1.0])
    plt.title(Drug_list[i])

plt.subplots_adjust(wspace=0.5, hspace=0.8) 
plt.show()	
