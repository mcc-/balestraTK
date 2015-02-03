balestraTK
==========

Balestra ToolKit

Parser for drug-target interaction related databases


DrugBank parser usage:
from balestratk.drugbank import *
db = DrugBank('/mnt/data/drugbank_v4.1/') #Folder where you keep DrugBank data
print db['aspirin'].synonyms #should work.
