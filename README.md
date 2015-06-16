BalestraTK: Balestra ToolKit
==========

Parser for drug-target interaction related databases.
There are two main uses of BalestraTK: 

1. To parse DrugBank
2. To parse STITCH-UniProt-PubChem data

If the user prefers a small but more richly annotated set of interactions,
drugs and targets then DrugBank should be preferred. If the user prefers a
large number of interactions then the STITCH-based interactions dataset should
be preferred. 

The stats for the two databases are as follows:

| Dataset | No of Drugs    | No of Proteins | No of Interactions |
| --------|---------|-------|------|
| DrugBank (v4.1)  | 7,740   | 4,103    | 15,120 |
| STITCH (v4.0; human) | 141,799 | 19,488 | 4,523,609 |

To use the toolkit, the users must first make sure they have the following
dependencies: 
- GraphLab Create: https://dato.com/products/create/
- Levenshtein: https://pypi.python.org/pypi/python-Levenshtein/

The usage logic is as follows: the user decides on a folder to keep the data,
the toolkit builds the data structures it needs to run in the first run from
that folder, in subsequent runs when the user points to that folder the code
automatically retrieves the previously built data structures and initializes
rapidly. 

The very first use from a folder triggers a necessary building step where
the data is parsed for easier access in later uses. From the user's point of
view, this means that there is a one-time time cost only on the very first use,
which enables very rapid use at all the following uses.

To use the toolkit for DrugBank the user must first have downloaded DrugBank, 
specifically, `drugbank.xml` available from http://www.drugbank.ca/downloads 
with v4.1 being the latest version available for download at the time of writing 
thus it is the most recent version with which the provided parser is guaranteed
to work with. The user is required to download the data because DrugBank license 
explicitly states that the redistribution of DrugBank v4.1 data is not allowed. 
After downloading, the data of DrugBank v4.1 can be processed by BalestraTK 
simply by pointing to the folder with the data. 

```python
from balestratk.drugbank import *
db = DrugBank('/folder/with/drugbank/data/') #Folder where you keep DrugBank data
print db['aspirin'].synonyms # prints synonyms of aspirin in DrugBank
print d.targets # prints the targets of aspirin
print d.transporters # prints the transporters of aspirin

for d in db.iterDrugs(): 
    print d.name # iterates over the drugs

for p in db.iterProts():
    print p.name # iterates over the targets

# The following selects all approved drugs
approved_drugs = []
for d in db.iterDrugs():
    if 'approved' in d.groups:
        approved_drugs.append(d)

# The following iterates over all interactions of approved drugs
for d in approved_drugs:
    if getAttr(d, 'targets', None): #necessary as some drugs have no targets
                                    #returns targets if they exist, None o.w.
        print d.targets 
```

To use the toolkit to access STITCH data, the user must simply point at the
folder where they would like to keep the STITCH data. If this folder already
has the STITCH files, the code uses those. If not, then the necessary files are
automatically retrieved. This is because STITCH license, at the time of
writing, is permissive to this type of access. The usage style is exactly as above:

```python
from balestratk.stitch import *
st = STITCH('/folder/with/stitch/data/') #Folder where you keep STITCH data
c = st['aspirin'] # automatically retrieves corresponding record in O(1) time
print c.smiles # prints the structure of the chemical in SMILES form
for intr in st.iter_chem_intr(c):
    print intr # prints the interactions of aspirin

p = st['adrb2'] # automatically retrieves corresponding record in O(1) time
print p.sequence # prints the sequence of the protein
for intr in st.iter_prot_intr(p):
    print intr # prints the interactions of ADRB2

# Likewise all the interactions in STITCH can be iterated over as follows
for intr in st.iter_intr():
    print intr
```
