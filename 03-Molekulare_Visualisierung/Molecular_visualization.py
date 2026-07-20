#!/usr/bin/env python
# coding: utf-8

# # Visualisierung molekularer Daten und Strukturen mit **RDKit** und **py3Dmol**

# Die interaktiven 3D-Strukturen von Py3dMol sind nur beim Ausführen des Notebooks in Browser verfügbar. Auf GitHub werden sie als leere Felder angezeigt.

# #### **RDKit** importieren und installieren
# 
# RDKit in Python:
#     https://www.rdkit.org/docs/GettingStartedInPython.html
# 
# RDKit Cookbook:
#     https://www.rdkit.org/docs/Cookbook.html

# In[1]:


from rdkit import Chem
from rdkit.Chem import AllChem


# #### **py3Dmol** installieren und importieren:
# 
# Viewer: py3Dmol 
#     https://pypi.org/project/py3Dmol/
#     
# 
# py3Dmol Tuturial:
#     https://colab.research.google.com/drive/1T2zR59TXyWRcNxRgOAiqVPJWhep83NV_?usp=sharing#scrollTo=0N6PEYCyGx5g
# 
# 

# In[2]:


import py3Dmol


# ## RDKit: Darstellungen und Koordinaten von Molekülen aus Identifikatoren generieren

# #### Darstellungen von Molekülen aus SMILES und InChI generieren

# In[3]:


# Generate Molecule from SMILES: Example Toluene
ex1 = Chem.MolFromSmiles('Cc1ccccc1')
ex1


# In[4]:


# generate Example 1 (Toluene) from InChI
ex1_i = Chem.MolFromInchi('InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3')
ex1_i


# In[5]:


# Generate Molecule from SMILES: Example aspirin
ex2 = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
ex2


# In[6]:


# Generate Molecule from SMILES example Benzoic acid
ex3 = Chem.MolFromSmiles('C1=CC=C(C=C1)C(=O)O')
ex3


# #### Eigenschaften der Moleküle bekommen

# In[7]:


# Number of Atoms (Example 3)
n_at = ex3.GetNumAtoms()
print('number of atoms = ', n_at)


# In[8]:


# Number of Bonds (Example 3)
n_bonds = ex3.GetNumBonds()
print('number of bonds = ', n_bonds)


# In[9]:


# Generate InCHi from molecule 
ex3_inchi = Chem.inchi.MolToInchi(ex3)
print('InChi of molecule = ', ex3_inchi)


# In[10]:


# Regenerate Example 3 (Aspirin) from InChI as a test
ex3_i = Chem.MolFromInchi(ex3_inchi)
ex3_i


# #### 3D-Koordinaten eines Moleküls generieren und visualisieren

# In[11]:


# Generate coordinates and connectivity
print(Chem.MolToMolBlock(ex3))    


# In[12]:


# Add Hydrogens
ex3_H = Chem.AddHs(ex3)
ex3_H


# In[13]:


# Generate 3D Coordinates with hydrogens
params = AllChem.ETKDGv3()
params.randomSeed = 0xf00d # optional random seed for reproducibility
AllChem.EmbedMolecule(ex3_H, params)
ex3H_crd = Chem.MolToMolBlock(ex3_H)
print(ex3H_crd)  


# In[14]:


# Visualize with py3Dmol
xyzview = py3Dmol.view(width=400,height=400)
xyzview.addModel(ex3H_crd)
xyzview.setStyle({'stick':{}})
xyzview.zoomTo()
xyzview.show()


# In[15]:


# Generate Coordinates in PDB format
print(Chem.MolToPDBBlock(ex3_H,flavor=4))


# In[16]:


# Generate Coordinates in xzy format
print(Chem.MolToXYZBlock(ex3_H))


# In[17]:


# Store Coordinates in PDB format
Chem.MolToPDBFile(ex3_H, 'test.pdb', flavor= 4)


# In[18]:


# Store Coordinates in xyz format
Chem.MolToXYZFile(ex3_H, 'test.xyz')


# ## Load Molecules from File in py3Dmol

# Example coordinates of carbon nanotube:
#     https://www.ks.uiuc.edu/Training/Tutorials/science/nanotubes/files/
# 

# In[19]:


# First we load the xzy file we have generated before
with open("test.xyz") as ifile:
    sys = "".join([x for x in ifile])


# In[20]:


xyzview = py3Dmol.view(width=400,height=400)
xyzview.addModel(sys)
xyzview.setStyle({'stick':{}})
xyzview.zoomTo()
xyzview.show()


# #### Carbon Nanotubes aus externen Files laden
# 
# File Download: 
#     https://www.ks.uiuc.edu/Training/Tutorials/science/nanotubes/files/

# In[21]:


with open("nanotubes.pdb") as ifile:
    system = "".join([x for x in ifile])


# In[22]:


view = py3Dmol.view(width=600, height=600)
view.addModelsAsFrames(system)
view.setStyle({'model': -1}, {"stick": {'color': 'spectrum'}})
view.zoomTo()
view.show()


# In[23]:


view = py3Dmol.view(width=600, height=600)
view.addModelsAsFrames(system)
view.setStyle({'model': -1}, {"sphere": {'color': 'spectrum'}})
view.zoomTo()
view.show()


# In[24]:


view = py3Dmol.view(width=600, height=600)
view.addModelsAsFrames(system)
view.setStyle({'model': -1}, {"stick": {'color': 'spectrum'}})
view.addSurface(py3Dmol.VDW,{'opacity':0.7,'color':'spectrum'})
view.zoomTo()
view.show()


# In[25]:


view = py3Dmol.view(width=600, height=600)
view.setBackgroundColor(0x000000)
view.addModelsAsFrames(system)
view.setStyle({'model': -1}, {"stick": {'color': 'spectrum'}})
view.zoomTo()
view.show()


# ### Biomoleküle aus der from Protein Database visualisieren
# 
# Protein Database:
#         https://www.rcsb.org/
# 
# Beispiel: PDB Struktur 2OWM, Motor Protein:
#     https://www.rcsb.org/structure/2OWM

# In[26]:


view = py3Dmol.view(query='pdb:2owm')
view.setStyle({'cartoon': {'color':'spectrum'}})
view.show()


# In[27]:


view = py3Dmol.view(query='pdb:2owm')
chA = {'chain':'A'}
chB = {'chain':'B'}
chC = {'chain':'C'}
chD = {'chain':'D'}
view.setStyle(chA,{'cartoon': {'color':'white'}})
view.setStyle(chB,{'cartoon': {'color':'black'}})
view.setStyle(chC,{'cartoon': {'color':'yellow'}})
view.setStyle(chD,{'cartoon': {'color':'blue'}})
view.show()


# In[ ]:




