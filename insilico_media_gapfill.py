#!/usr/bin/env python
# coding: utf-8

# In[1]:


from cobra import Model, Reaction, Metabolite
from riptide import *
import pandas as pd
from copy import *
from cobra.medium import minimal_medium
from difflib import *


# In[4]:


# Protein sequence location https://www.uniprot.org/proteomes/UP000035707
model = cobra.io.read_sbml_model("gvaginalis_gapseq_draftmodel.xml")


# In[12]:


model_metabolite_names = []
for metabolite in model.metabolites:
    model_metabolite_names.append(metabolite.name)   


# In[13]:


#synthetic vaginal media 
vaginal_media_names = ["sodium",#NaCl originally
                        "potassium", #KCl originally
                        #K2HPO4 buffering agent not included
                        #KH2PO4 buffering agent not included
                        "Dextrin", #dextrose originally
                        "Cysteine", #cysteine HCl originally
                        "Glycogen", #"mucin",
                        #Tween 20 not included
                        "Urea", #"phytomenadione",
                        "Heme", #originally hemin
                        #"Albumin",
                        "Mg", #originally MgSO4
                        "sulfate", #originally MgSO4
                        "H2CO3", #originally NaHCO3
                        "Biotin",
                        "1-Phosphatidyl-myo-inositol", #originally myo-inositol
                        "niacin", #originally niacinamide
                        "pyridoxal", #pyridoxine HCl originally
                        "Thiamin",#thiamine HCl originally
                        "pantothenic acid", #originally D-Calcium pantothenate
                        "Ca2", #originally D-Calcium pantothenate
                        "Folic acid",
                        "4-aminobenzoate", #p-Aminobenzoic acid originally
                        "Choline", #choline chloride originally
                        "Riboflavin",
                        "L-ascorbic acid",
                        #"retinol",
                        "ergosterol", #calciferol
                        #"Cyanocobalamin",
                        "water"]
#KSFM --> based on BHI enriched media conditions from Deb
KSFM_media_names = ["alanine",
                    "aspartic acid",
                    "glutamic acid",
                    "histidine",
                    "leucine",
                    "methionine",
                    "proline",
                    "threonine",
                    "arginine",
                    "tyrosine",
                    "cysteine",
                    "glycine",
                    "isoleucine",
                    "lysine",
                    "phenylalanine",
                    "serine",
                    "tryptophan",
                    "valine",
                    "glucose",
                    "riboflavin",
                    "calcium pantothentate",
                    "folate",
                    "niacin",
                    "pyridoxine",
                    "biotin",
                    "cobalt",
                    "sodium",
                    "chloride",
                    "potassium",
                    "phosphate",
                    "calcium",
                    "magnesium",
                    "Fe3+",
                    "manganese",
                    "carbonate",
                    "H2O",
                    "HYXN", #originally hypoxanthine
                    "thyminose",
                    "cytosine","uracil"]
#NYCIII media --> ATCC medium: 1685 NYC III
NYCIII_media_names = ["sodium",
                      "chloride",
                      "H2O",
                      "D-glucose",
                      #hepes buffering agent not included
                      #fresh yeast extract components:
                      #https://pubs.acs.org/doi/full/10.1021/acsfoodscitech.0c00131?ref=recommended&
                      "biotin",
                      "thiamin",
                      "riboflavin",
                      "niacin",
                      "pantothenic acid",
                      "glutamate",
                      "glycine",
                      "alanine",
                      "valine",
                      "folic acid",
                      "potassium",
                      "zinc",
                      "Fe3+",
                      #"glutathione",
                      #"selenium",
                      "glucose", #originally glucan
                      #proteose peptone --> primarily source of amino acids & nitrogen
                      #"nitrogen",
                      "arginine",
                      "asparagine",
                      "aspartic acid",
                      "cysteine",
                      "glutamine",
                      "glutamic acid",
                      "glycine",
                      "histidine",
                      "isoleucine",
                      "leucine",
                      "lysine",
                      "methionine",
                      "phenylalanine",
                      "proline",
                      "serine",
                      "threonine",
                      "tryptophan",
                      "tyrosine",
                      "valine",
                      #horse serum composition: https://www.mdpi.com/2218-1989/10/7/298/htm
                      #"3-Hydroxybutyrate",
                      #"2-Hydroxyisobutyrate",
                      "Acetoacetate",
                      "acetone",
                      "citrate",
                      "creatine",
                      #"diemthyl sulfone",
                      "diemthylglycine",
                      "fumarate",
                      "glutamine",
                      "glycerol",
                      "lactate",
                      "mannose",
                      "methanol",
                      "myo-inositol",
                      "pyruvate",
                      "succinate",
                      "trimethylamine",
                      "formate",
                      "arabinose",
                      "betaine",
                      #"creatinine",
                      #"sarcosine",
                      #"diemthylamine",
                      "acetate",
                      "2-hydroxybutyrate",
                      "ethanol",
                      "3-hydroxyisobutyrate"]


# In[33]:


#now that we have the list of exchange reaction IDs we can go ahead and take care of our in silico media 
#and gapfill so that our model can grow on our media
SVM_metabolite_id = ["cpd00971",
                        "cpd00205",
                        "cpd11594",
                        "cpd00084",
                        "cpd00155",
                        "cpd00073",
                        "cpd00254",
                        "cpd00048",
                        "cpd00242",
                        "cpd00028",
                        "cpd00104",
                        "cpd11822",
                        "cpd00218",
                        "cpd00215",
                        "cpd00305",
                        "cpd00644",
                        "cpd00393",
                        "cpd00443",
                        "cpd00098",
                        "cpd00220",
                        "cpd00059",
                        "cpd01170",
                        "cpd00001",
                        "cpd00063"]
KSFM_metabolite_id =["cpd00035",
                        "cpd00041",
                        "cpd00023",
                        "cpd00119",
                        "cpd15143",
                        "cpd00060",
                        "cpd00129",
                        "cpd00161",
                        "cpd00051",
                        "cpd00069",
                        "cpd00547",
                        "cpd00033",
                        "cpd00322",
                        "cpd00039",
                        "cpd00066",
                        "cpd00054",
                        "cpd00065",
                        "cpd00156",
                        "cpd00027",
                        "cpd00220",
                        "cpd00644",
                        "cpd00393",
                        "cpd00218",
                        "cpd00215",
                        "cpd00104",
                        "cpd00149",
                        "cpd00971",
                        "cpd00099",
                        "cpd00205",
                        "cpd00009",
                        "cpd00063",
                        "cpd00254",
                        "cpd10516",
                        "cpd00030",
                        "cpd00242",
                        "cpd00001",
                        "cpd00226",
                        "cpd01242",
                        "cpd00307",
                        "cpd00092"]
NYCIII_metabolite_id = ["cpd00971",
                        "cpd00099",
                        "cpd00001",
                        "cpd00027",
                        "cpd00104",
                        "cpd00305",
                        "cpd00220",
                        "cpd00218",
                        "cpd00644",
                        "cpd00393",
                        "cpd00023",
                        "cpd00033",
                        "cpd00035",
                        "cpd00156",
                        "cpd00205",
                        "cpd10516",
                        "cpd00034",
                        "cpd00027",
                        "cpd00051",
                        "cpd00132",
                        "cpd00041",
                        "cpd00547",
                        "cpd00053",
                        "cpd00023",
                        "cpd00033",
                        "cpd00119",
                        "cpd00322",
                        "cpd00039",
                        "cpd15143",
                        "cpd00060",
                        "cpd00066",
                        "cpd00054",
                        "cpd00065",
                        "cpd00129",
                        "cpd00161",
                        "cpd00069",
                        "cpd00156",
                        "cpd03561",
                        "cpd00876",
                        "cpd00142",
                        "cpd00178",
                        "cpd00137",
                        "cpd00281",
                        "cpd00106",
                        "cpd00253",
                        "cpd00100",
                        "cpd00159",
                        "cpd00138",
                        "cpd00116",
                        "cpd00121",
                        "cpd00020",
                        "cpd00036",
                        "cpd00441",
                        "cpd00047",
                        "cpd00224",
                        "cpd00540",
                        "cpd01550",
                        "cpd00029",
                        "cpd00363"]
#create list of synthetic vaginal media exchange reactions
SVM_open_exchanges = []
for metabolite in SVM_metabolite_id:
    for exchange in model.exchanges:
        if metabolite in exchange.reaction:
            SVM_open_exchanges.append(exchange.id)

#create list of KSFM exchange reactions
KSFM_open_exchanges = []
for metabolite in KSFM_metabolite_id:
    for exchange in model.exchanges:
        if metabolite in exchange.reaction:
            KSFM_open_exchanges.append(exchange.id)

#create list of NYC III exchange reactions
NYCIII_open_exchanges = []
for metabolite in NYCIII_metabolite_id:
    for exchange in model.exchanges:
        if metabolite in exchange.reaction:
            NYCIII_open_exchanges.append(exchange.id)


# In[ ]:


#we used gapseq to gapfill our model so that it can grow on all three media conditions
model_gapfilled = cobra.io.read_sbml_model("gvaginalis_gapfilled.xml")


# In[ ]:


#Create list of minimal media exchange reactions
min_media = minimal_medium(model_gapfilled, 0.5, open_exchanges=True)
minmedia_open_exchanges = list(min_media.index)


# In[16]:


def changeMedia(model, media, limEX=[]):
    modelOutput = deepcopy(model_gapfilled)

    # Set the new media conditions
    for ex in modelOutput.exchanges:
        ex.upper_bound = 1000
        ex.lower_bound = 0
        
    # Minimal Media
    if media == 1:
        for exchange in modelOutput.reactions:
            if exchange.id in minmedia_open_exchanges:
                exchange.lower_bound = -1000
        
    # Synthetic Vaginal Media
    elif media == 2:
        for exchange in modelOutput.reactions:
            if exchange.id in SVM_open_exchanges:
                exchange.lower_bound = -1000
        
    # Minimal Media + Synthetic Vaginal Media 
    elif media == 3:
        for exchange in modelOutput.reactions:
            if exchange.id in minmedia_open_exchanges:
                exchange.lower_bound = -1000
            elif exchange.id in SVM_open_exchanges:
                exchange.lower_bound = -1000
                
     # NYCIII
    elif media == 4:
        for exchange in modelOutput.reactions:
            if exchange.id in NYCIII_open_exchanges:
                exchange.lower_bound = -1000   
    # KSFM
    elif media == 5:
        for exchange in modelOutput.reactions:
            if exchange.id in KSFM_open_exchanges:
                exchange.lower_bound = -1000
   
    elif media == 0:
        print('all exchange bounds set to [0,1000]')
            
    else:
        print('unrecognized media condition. Please enter 1 for Minimal Media; 2 for Synthetic Vaginal Media; 3 for minimal media; 4 for NYCIII; 5 for KSFM')
   
    return(modelOutput) 

