#!/usr/bin/env python
# coding: utf-8

# In[1]:


from cobra import Model, Reaction, Metabolite
from riptide import *
import pandas as pd
from copy import *
from cobra.medium import minimal_medium
from difflib import *


# In[2]:


# Protein sequence location https://www.uniprot.org/proteomes/UP000035707
model = cobra.io.read_sbml_model("gvaginalis.xml")


# In[3]:


solution = model.optimize()
print(solution)


# In[4]:


model_metabolite_names = []
for metabolite in model.metabolites:
    model_metabolite_names.append(metabolite.name)   


# In[5]:


#synthetic vaginal media 
vaginal_media_names = ["sodium",#NaCl originally
                        "potassium", #KCl originally
                        #K2HPO4 buffering agent not included
                        #KH2PO4 buffering agent not included
                        "Dextrose",
                        "Cysteine", #cysteine HCl originally
                        "Glycogen", "mucin",
                        #Tween 20 not included
                        "Urea", "phytomenadione", "Hemin", "Albumin",
                        "magnesium", #originally MgSO4
                        "sulfate", #originally MgSO4
                        "bicarbonate", #originally NaHCO3
                        "Biotin", "myo-Inositol",
                        "niacinamide",
                        "Pyridoxine", #pyridoxine HCl originally
                        "Thiamine",#thiamine HCl originally
                        "D-Calcium pantothenate", "Folic acid", "p-Aminobenzoic acid",
                        "Choline chloride", "Riboflavin", "L-ascorbic acid", "retinol", "calciferol",
                        "Cyanocobalamin","water"]
metabolite_match = []
for metabolite in vaginal_media_names:
    metabolite_match.append(get_close_matches(metabolite, model_metabolite_names))

metabolite_match = set([item for sublist in metabolite_match for item in sublist])

#review matches to remove any incorrectly matched metabolites

SVM_metabolite_id = []
for metabolite in model.metabolites:
    if metabolite.name in metabolite_match:
        SVM_metabolite_id.append(metabolite.id)


# In[6]:


#KSFM --> based on BHI enriched media conditions from Deb
KSFM_media_names = ["alanine","aspartic acid","glutamic acid","histidine","leucine","methionine","proline",
                    "threonine","arginine","tyrosine","cysteine","glycine","isoleucine","lysine","phenylalanine",
                    "serine","tryptophan","valine","glucose","riboflavin","calcium pantothentate","folate",
                    "niacin","pyridoxine","biotin","cobalt","sodium","chloride","potassium","phosphate","calcium",
                    "magnesium","Iron (Fe3+)","manganese","carbonate","H2O H2O","hypoxanthine","thyminose","cytosine","uracil"]
metabolite_match = []
for metabolite in KSFM_media_names:
    metabolite_match.append(get_close_matches(metabolite, model_metabolite_names))

metabolite_match = set([item for sublist in metabolite_match for item in sublist])

#review matches to remove any incorrectly matched metabolites

KSFM_metabolite_id = []
for metabolite in model.metabolites:
    if metabolite.name in metabolite_match:
        KSFM_metabolite_id.append(metabolite.id)


# In[7]:


#NYCIII media --> ATCC medium: 1685 NYC III
NYCIII_media_names = ["sodium","chloride", "H2O H2O", "glucose",
                        #hepes buffering agent not included
                        #fresh yeast extract components:
                        #https://pubs.acs.org/doi/full/10.1021/acsfoodscitech.0c00131?ref=recommended&
                        "biotin", "thiamine", "riboflavin", "niacin","pantothenic acid", "glutamic acid","glycine",
                        "alanine","valine","folic acid", "potassium", "zinc", "iron (Fe3+)", "sodium", "glutathione",
                        "selenium", "glucan",
                        #proteose peptone --> primarily source of amino acids & nitrogen
                        "nitrogen", "arginine", "asparagine", "aspartic acid", "cysteine", "glutamine", "glutamic acid",
                        "glycine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline",
                        "serine", "threonine", "tryptophan", "tyrosine", "valine",
                        #horse serum composition: https://www.mdpi.com/2218-1989/10/7/298/htm
                        "3-Hydroxybutyrate", "2-Hydroxyisobutyrate", "Acetoacetate", "acetone","citrate", "creatine",
                        "diemthyl sulfone", "diemthylglycine", "fumarate", "glutamine", "glycerol", "lactate",
                        "mannose", "methanol", "myo-inositol", "pyruvate", "succinate", "trimethylamine", "formate",
                        "phenylalanine", "arabinose", "betaine", "creatinine", "sarcosine", "diemthylamine", "acetate",
                        "2-hydroxybutyrate", "ethanol", "3-hydroxyisobutyrate"]
metabolite_match = []
for metabolite in NYCIII_media_names:
    metabolite_match.append(get_close_matches(metabolite, model_metabolite_names))

metabolite_match = set([item for sublist in metabolite_match for item in sublist])

NYCIII_metabolite_id = []
for metabolite in model.metabolites:
    if metabolite.name in metabolite_match:
        NYCIII_metabolite_id.append(metabolite.id)


# In[8]:


#create list of existing exchange ids
media_ids = list(NYCIII_metabolite_id + KSFM_metabolite_id + SVM_metabolite_id)
media_ids = list(set([x[:-2] for x in media_ids]))
media_exchange_ids = ["EX_"+ x + "_e" for x in media_ids]
model_exchange_ids = []
for exchange in model.exchanges:
    if exchange.id in media_exchange_ids:
        model_exchange_ids.append(exchange.id)


# In[10]:


#Find missing exchange reactions 
exchange_difference = [x for x in media_exchange_ids if x not in model_exchange_ids]
missing_exchange_metID = [x[3:] for x in exchange_difference]

#add missing extracellular metabolites
for cpd in missing_exchange_metID:
    metabolite = Metabolite(cpd)
    for ids in model.metabolites:
        if cpd in ids.id:
            metabolite.name = ids.name
    metabolite.compartment = "extracellular"
    model.add_metabolites([metabolite])


# In[11]:


#add missing exchange reactions
for metabolite in missing_exchange_metID:
    reaction = Reaction('EX_' + metabolite)
    reaction.name = model.metabolites.get_by_id(metabolite).name + 'exchange'
    reaction.subsystem = 'exchange'
    reaction.lower_bound = 0 
    reaction.upper_bound = 1000
    reaction.add_metabolites({model.metabolites.get_by_id(metabolite):-1.0})
    model.add_reactions([reaction])


# In[12]:


#identify what metabolites have transport reactions
#identify which reactions have the stripped media ids >= 2 times in their reaction formula
transport_ids = []
for reaction in model.reactions: 
    for formula in reaction.reaction:
        for ids in media_ids:
            if formula.count(ids) >= 2:
                transport_ids.append(reaction.id)


# In[ ]:


#To Do
#4. add missing transport reactions 


# In[ ]:


#now that we have the list of exchange reaction IDs we can go ahead and take care of our in silico media 
#and create a function to change our media conditions


# In[13]:


#Create list of minimal media exchange reactions
min_media = minimal_medium(model, 0.5, open_exchanges=True)
minmedia_open_exchanges = list(min_media.index)

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


# In[17]:


def changeMedia(model, media, limEX=[]):
    modelOutput = deepcopy(model)

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


# In[ ]:


practice = changeMedia(model,1)

