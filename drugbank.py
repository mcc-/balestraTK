from xml.etree import ElementTree as ET

import shelve
import zipfile
import os

def clrns(text):
    return text[text.find('}')+1:].replace('-','')

def combine(l):
    if set([type(p) for p in l]).issubset({dict, type(None)}):
        r = {}
        for e in l:
            if type(e) == dict:
                for k,v in e.iteritems():
                    if k in r:
                        if type(r[k])==list:
                            r[k].append(v)
                        else:
                            r[k] = [r[k], v]
                    else:
                        r[k] = v
    elif set([type(p) for p in l]).issubset({dict, type(None), list}):
        r = []; rd = []
        for e in l:
            if type(e) == list:
                r.append(combine(e))
            elif type(e) == dict:
                rd.append(e)
        if len(rd) > 0:
            r.append(combine(rd))
    else:
        if len(l) == 1:
            return l[0]
        else:
            return l
    return r

def getMapFromElement(elmt,parent=None):
    children = elmt.getchildren()
    tag = clrns(elmt.tag)
    if len(children) == 0:
        if elmt.text is not None:
            if (parent is not None)  and  (tag+'s' == parent):
                return elmt.text.strip()
            else:
                return {tag:elmt.text.strip()}
        else:
            return {tag:None}
    else:
        res = combine([getMapFromElement(child,parent=tag) for child in children])
        return res

class Drug():
    def __init__(self,elmt):
        for child in elmt.getchildren():
            tag = clrns(child.tag)
            grandchildren = child.getchildren()
            if len(grandchildren) == 0 and child.text is not None:
                if getattr(self,tag,None) == None or 'primary' in child.attrib:
                    setattr(self, tag, child.text.strip())
            elif len(grandchildren) > 0:
                if tag == 'atccodes':
                    self.atc = grandchildren[0].attrib['code']
                elif tag == 'groups':
                    setattr(self,tag,[g.text for g in grandchildren])
                elif tag in ['targets', 'transporters', 'enzymes']:
                    setattr(self,tag,[Protein(g) for g in grandchildren])
                elif tag == 'categories':
                    self.categories = [g.getchildren()[0].text for g in grandchildren]
                else:
                    setattr(self, tag, getMapFromElement(child))
        self.identifier = self.drugbankid
        del(self.drugbankid)

    def __repr__(self):
        return '<Drug: '+self.identifier+' '+self.name+'>'

class Protein():
    def __init__(self,elmt):
        for child in elmt.getchildren():
            tag = clrns(child.tag)
            grandchildren = child.getchildren()
            if len(grandchildren) == 0 and child.text is not None:
                if getattr(self,tag,None) == None or 'primary' in child.attrib:
                    setattr(self, tag, child.text.strip())
            elif len(grandchildren) > 0:
                if tag == 'polypeptide':
                    for gc in grandchildren:
                        gctag = clrns(gc.tag)
                        element = getMapFromElement(gc)
                        if type(element)==dict and                             len(element)==1 and                             element.keys()[0]==gctag:
                            setattr(self, gctag, element.values()[0])
                        else:
                            setattr(self, gctag, element)
                else:
                    setattr(self, tag, getMapFromElement(child))
        self.identifier = self.id
        del(self.id)

    def __repr__(self):
        return '<Protein:'+self.identifier+' '+self.name+'>'

class DrugBank():
    def __init__(self,folder,rebuild=False):
        "*folder* is the folder with DrugBank v4 XML file to parse"
        if rebuild is True or not(os.path.isfile(folder+os.sep+'drugs.shl')
                              and os.path.isfile(folder+os.sep+'prots.shl')):
            xmlfile = zipfile.ZipFile(folder+os.sep+'drugbank.xml.zip').open('drugbank.xml')
            xml = ET.parse(xmlfile)
            root = xml.getroot()
            self.drugs = shelve.open(folder+os.sep+'drugs.shl',protocol=2,flag='n',writeback=True)
            self.prots = shelve.open(folder+os.sep+'prots.shl',protocol=2,flag='n',writeback=True)
            self.actions = shelve.open(folder+os.sep+'actions.shl',protocol=2,flag='n',writeback=True)
            self.info = shelve.open(folder+os.sep+'info.shl',protocol=2,flag='n',writeback=True)
            self.info['no_drug'] = 0
            self.info['no_prot'] = 0
            self.info['no_intr'] = 0
            self.name_map = {}
            for child in root.getchildren():
                self.__registerDrug(child)
            self.info['name_map'] = self.name_map
            self.drugs.sync()
            self.prots.sync()
            self.actions.sync()
            self.info.sync()
        else:
            self.drugs = shelve.open(folder+os.sep+'drugs.shl',protocol=2,flag='r')
            self.prots = shelve.open(folder+os.sep+'prots.shl',protocol=2,flag='r')
            self.actions = shelve.open(folder+os.sep+'actions.shl',protocol=2,flag='r')
            self.info = shelve.open(folder+os.sep+'info.shl',protocol=2,flag='r')
            self.name_map = self.info['name_map']


    def __registerDrug(self,elmt):
        "*elmt* is the ElementTree instance of drug"
        d = Drug(elmt)
        names = [d.name]; name_sets = []
        if getattr(d,'synonyms',None) is not None: name_sets.append(d.synonyms)
        if getattr(d,'brands',None) is not None: name_sets.append(d.brands)
        for name_set in name_sets:
            if type(name_set == list):
                names.extend(name_set)
            elif type(name_set == str):
                names.append(name_set)
        for name in names:
            self.name_map[name.upper()] = d.identifier
        self.info['no_drug'] = self.info['no_drug'] + 1
        actions = {}
        if getattr(d,'targets',None) is not None:
            self.info['no_intr'] = self.info['no_intr'] + len(d.targets)
            for t in d.targets:
                actions.update({t.identifier:t.knownaction})
                self.__registerProt(t)
        self.actions[d.identifier] = actions
        self.drugs[d.identifier] = d


    def __registerProt(self,prot):
        "*prot* is the ElementTree instance of protein"
        if prot.identifier not in self.prots:
            self.info['no_prot'] = self.info['no_prot'] + 1
            names = [prot.name]
            genename = getattr(prot,'genename',None)
            if genename is not None: names.append(genename)
            if getattr(prot,'synonyms',None) is not None:
                if type(prot.synonyms == list):
                    names.extend(prot.synonyms)
                elif type(prot.synonyms == str):
                    names.append(prot.synonyms)
            if getattr(prot,'externalidentifiers',None) is not None:
                if 'identifier' in prot.externalidentifiers:
                    ids = prot.externalidentifiers['identifier']
                    if type(ids) == list:
                        names.extend(ids)
                    else:
                        names.append(ids)
            for name in names:
                self.name_map[name.upper()] = prot.identifier
            del prot.knownaction
            self.prots[prot.identifier] = prot


    def __getitem__(self,key):
        """ Check for the name in the following order:
        Drug ID,
        Protein ID,
        Drug Name,
        Protein Name
        Then check if the input is a prefix of a name"""
        key=key.upper()
        try:
            return self.drugs[key]
        except KeyError:
            pass
        try:
            return self.prots[key]
        except KeyError:
            pass
        try:
            return self.__getitem__(self.name_map[key])
        except KeyError:
            pass
        res = []
        for k in self.name_map.iterkeys():
            if k.startswith(key):
                res.append(k)
        return res


    def iterDrugs(self):
        for d in self.drugs.itervalues():
            yield d

    def iterProts(self):
        for p in self.prots.itervalues():
            yield p

    def __repr__(self):
        return "<DrugBank {0:g} drugs, {1:g} proteins, {2:g} links>".format(                                                                    self.info['no_drug'],                                                                    self.info['no_prot'],                                                                    self.info['no_intr'])

