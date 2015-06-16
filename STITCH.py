import multiprocessing
import gzip
import re
import json
import shelve
import cPickle as pickle
import httplib, urllib2
from collections import defaultdict
import sys,os
import graphlab as gl
import graphlab.aggregate as agg
import numpy as np
import Levenshtein

class Protein:
    def __init__(self, stitchid=None):
        if stitchid is None:
            raise ValueError('Must be instantiated with UniProt identifier')
        self.acc = []
        self.godict = defaultdict(list)
        self.goterms = []
        self.sequence = ''
        try:
            self.organism, protid = stitchid.split('.')
        except ValueError:
            protid = stitchid
            self.organism = "9606"
        self.stitchid=self.organism+'.'+protid
        self.name = protid
        self.fullname = protid
        self.pfam = []
        self.uniprotid = protid
        try:
            tickets = 5; success = False
            while not(success) and tickets > 0:
                try:
                    res = urllib2.urlopen("http://www.uniprot.org/uniprot/?query=taxonomy%3a"+self.organism+"+AND+"\
                        +protid+"&format=txt")
                    success = True
                except httplib.BadStatusLine:
                    tickets -= 1
                    time.sleep(0.5)
            if success:
                name_pattern = re.compile(r"[ ]?(?P<categ>\w+)=(?P<val>(\w+)(, (\w+))*)")
                seq_mode = False
                for line_raw in res:
                    line = [s.strip(';') for s in line_raw.strip().split()]
                    if line[0] == "ID":
                        self.uniprotid = line[1]
                    if line[0] == "AC":
                        self.acc.extend(line[1:])
                    if line[0] == "DE":
                        if line[1] == "RecName:":
                            self.fullname = " ".join(line[2:])
                            if self.fullname.startswith("Full="):
                                self.fullname = self.fullname[5:]
                    if line[0] == "GN":
                        sgmts = line_raw[5:].strip().split(';')
                        try:
                            for sgmt in sgmts:
                                match = name_pattern.search(sgmt)
                                if match:
                                    resd = match.groupdict()
                                    setattr(self,resd['categ'].lower(),resd['val'].split(', '))
                        except IndexError:
                            print 'IndexError', stitchid
                    if line[0] == "DR":
                        try:
                            if line[1] == "GO":
                                self.godict[line[3][0]].append((line[2],line[3][2:]))
                                self.goterms.append(" ".join(line[2:]))
                            if line[1] == "Pfam":
                                self.pfam.append(line[2])
                        except IndexError:
                            print 'IndexError', stitchid
                    if seq_mode == True:
                        if line_raw[:2] == '  ':
                            self.sequence += ''.join(line)
                        else:
                            seq_mode = False
                    if line[0] == "SQ":
                        seq_mode = True
                if type(self.name) is list and len(self.name)>0:
                    self.name = self.name[0]
                elif type(self.name) is not str:
                    self.name = protid
            else:
                raise ValueError
        except ValueError:
            print "Error parsing {0:s}".format(protid)

    def __repr__(self):
        return "<Protein {0:s}>".format(self.name)


class STITCH:
    def _get_frame(self,fname,url):
        if os.path.isdir(self.folder+fname+'.gl'):
            return gl.load_sframe(self.folder+fname+'.gl')
        else:
            if fname.endswith('.gz') and os.path.isfile(self.folder+fname[:-3]):
                frame = gl.SFrame.read_csv(self.folder+fname[:-3], delimiter='\t')
            elif os.path.isfile(self.folder+fname):
                frame = gl.SFrame.read_csv(self.folder+fname, delimiter='\t')
            else:
                urllib2.urlopen(url)
                print 'Downloading data from STITCH:', fname
                with file(self.folder+fname,'wb') as f:
                    f.write(urllib2.urlopen(url).read())
                frame = gl.SFrame.read_csv(self.folder+fname, delimiter='\t')
                os.remove(self.folder+fname)
            frame.save(self.folder+fname+'.gl')
            return frame

    def __init__(self, stitch_folder=None,rebuild=False):
        if stitch_folder is None:
            stitch_folder = os.getcwd()
        if stitch_folder[-1] != os.sep:
            stitch_folder += os.sep
        self.folder = stitch_folder
        # Get interactions
        self.st = self._get_frame("9606.protein_chemical.links.detailed.v4.0.tsv.gz",
                             "http://stitch.embl.de/download/protein_chemical.links.detailed.v4.0/9606.protein_chemical.links.detailed.v4.0.tsv.gz")
        # Get chemicals
        self.chem = self._get_frame("chemicals.v4.0.tsv.gz",
                             "http://stitch.embl.de/download/chemicals.v4.0.tsv.gz")
        # Get aliases
        self.chem_aliases = self._get_frame("chemicals.aliases.v4.0.tsv.gz",
                             "http://stitch.embl.de/download/chemical.aliases.v4.0.tsv.gz")
        if rebuild == False and os.path.isfile(self.folder+'chemicals.shl'):
            self.chem_map = shelve.open(self.folder+'chemicals.shl',flag='r',protocol=2)
            self.name_map = shelve.open(self.folder+'name_map.shl',flag='r',protocol=2)
            self.prot_map = shelve.open(self.folder+'proteins.shl',flag='r',protocol=2)
        else:
            print 'Built STITCH files not found in {0:s}\nBuilding STITCH...'.format(self.folder)
            uniq_chem = st.st['chemical'].unique()
            self.chem = self.chem.filter_by(uniq_chem, 'chemical')
            self.chem.save(self.folder+"chemicals.v4.0.tsv.gz.gl")
            self.chem_aliases = self.chem_aliases.filter_by(uniq_chem, 'chemical')
            self.chem_aliases.save(self.folder+"chemicals.aliases.v4.0.tsv.gz.gl")
            self.name_map = shelve.open(self.folder+'name_map.shl',flag='n',protocol=2)
            self._build_chem()
            self._build_prot()
        self.n_intr = self.st.num_rows()
        self.n_drug = len(self.chem_map.keys())
        self.n_prot = len(self.prot_map.keys())

    def _build_chem(self):
        self.chem_map = shelve.open(self.folder+'chemicals.shl',flag='n',protocol=2)
        print 'Integrating chemical info indexed on PubChem CID...'
        l = self.chem.num_rows()
        for i, row in enumerate(self.chem):
            self.name_map[row['name'].lower()] = row['chemical']
            self.chem_map[row['chemical']] = Chemical(row['chemical'],
                                                      row['name'],
                                                      row['SMILES_string'],
                                                      row['molecular_weight'])
            if i%1000==0:
                print 'Completed {0:d} / {1:d}'.format(i,l)
        for row in self.chem_aliases:
            self.name_map[row['alias'].lower()] = row['chemical']
        self.name_map.sync()
        self.chem_map.sync()

    def _build_prot(self):
        self.prot_map = shelve.open(self.folder+'proteins.shl',flag='n',protocol=2)
        uniq_prot = self.st['protein'].unique(); l = uniq_prot.size()
        print 'Integrating protein info indexed on Ensembl ID...'
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()*4)
        for ind, p in enumerate(pool.imap_unordered(Protein,uniq_prot,chunksize=20)):
            self.prot_map[p.stitchid] = p
            self.name_map[p.name.lower()] = p.stitchid
            if ind%100==0:
                print 'Completed {0:d} / {1:d}'.format(ind,l)
        self.name_map.sync()
        self.prot_map.sync()

    def iter_intr(self):
        "Returns iterator over all interactions"
        for intr in self.st:
            yield intr

    def iter_chem_intr(self,chem):
        "Returns the interactions of chemical *chem*"
        for intr in self.st.filter_by(chem.cid,'chemical'):
            yield intr

    def iter_prot_intr(self,prot):
        "Returns the interactions of protein *prot*"
        for intr in self.st.filter_by(prot.stitchid,'protein'):
            yield intr

    def iter_chem(self):
        "Returns an iterator over the chemicals"
        for chem in self.chem_map.values():
            yield chem

    def iter_prot(self):
        "Returns an iterator over the proteins"
        for prot in self.prot_map.values():
            yield prot

    def __getitem__(self,k):
        if type(k) is not str:
            print "Input must be string, instead received {0:s}".format(k)
            return None
        try:
            return self.chem_map[k.upper()]
        except KeyError:
            pass
        try:
            return self.prot_map[k.upper()]
        except KeyError:
            pass
        try:
            r = self.name_map[k.lower()]
            return self.__getitem__(r)
        except KeyError:
            pass
        try:
            r_list = []
            for n in self.name_map.iterkeys():
                if Levenshtein.ratio(n,k.lower()) > 0.95:
                    r_list.append(self.__getitem__(self.name_map[n]))
            return r_list
        except KeyError:
            pass
        print "Query not found"
        return None

    def __repr__(self):
        return "<STITCH {0:d} drugs, {1:d} proteins, {2:d} links>".format(self.n_drug, self.n_prot, self.n_intr)


class Chemical:

    def __init__(self, cid, name, smiles, molweight):
        self.cid = cid
        self.name = name
        self.smiles = smiles
        self.molweight = molweight

    def __repr__(self):
        return '<Chemical {0:s}:{1:s}>'.format(self.cid,self.name)
