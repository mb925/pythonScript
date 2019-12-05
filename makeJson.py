import urllib.request, json, os

unipIds = []

# getting a list of all uniprots ids
def getId():
    response = urllib.request.urlopen("http://repeatsdb.bio.unipd.it/ws/search?query=average_unit:1TO9999999999&collection=uniprot_protein&show=uniprotid")
    data = json.load(response)

    for elem in data:
        unipIds.append(elem['uniprotid'])

def getSeq(id):
    try:
        response = urllib.request.urlopen("https://www.ebi.ac.uk/proteins/api/proteins/" + id)
        data = json.load(response)

        seq = data['sequence']['sequence']
        return seq
    except:
        print ('missing data for : https://www.ebi.ac.uk/proteins/api/proteins/' + id +' ')
        return -1


def getPdbsData(id):
    try:

        response = urllib.request.urlopen("https://www.ebi.ac.uk/pdbe/api/mappings/" + id)
        data = json.load(response)
        return data
    except:
        print('missing data for https://www.ebi.ac.uk/pdbe/api/mappings/' + id + '')
    return -1
def getPdbs(id, seq):
        data = getPdbsData(id)
        if data != -1:
            pdbs = {}
            for elem in data[id]['PDB']:
               pdbs[elem] = {}
               resObj = getResObj(elem)
               if resObj == -1:
                   print('missing data for https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/'+ elem +'')
                   continue
               for ch in data[id]['PDB'][elem]:
                   chain = ch['chain_id']
                   eId = ch['entity_id']
                   pdbs[elem][chain] = {}
                   resNum = ch['start']['residue_number']
                   shift = ch['unp_start'] - resNum
                   pdbs[elem][chain]['unp_start'] = ch['unp_start']
                   pdbs[elem][chain]['unp_end'] = ch['unp_end']
                   residues = getResList(resObj, elem, chain, eId)
                   if residues == None:
                       print('missing residues, skipping: ' + elem + chain + str(eId) +'')
                       continue
                   unpAut = unpToAut(seq, residues, shift)
                   pdbs[elem][chain]['pdb_indexes'] = unpAut
                   pdbs[elem][chain]['repeats_db'] = {}
                   entities = getPdbsEntities(elem, chain)
                   autUnp = autToUnp(ch['unp_start'], ch['unp_end'], residues, shift)
                   units = elemToUnp(entities[0], autUnp)
                   insertions = elemToUnp(entities[1], autUnp)
                   pdbs[elem][chain]['repeats_db']['units'] = units
                   pdbs[elem][chain]['repeats_db']['insertions'] = insertions

                   if len(entities[0]) != 0:
                        pdbs[elem][chain]['repeats_db']['start'] = entities[0][0][0]
                        pdbs[elem][chain]['repeats_db']['end'] = entities[0][-1][1]
                        pdbs[elem][chain]['repeats_db']['classification'] = entities[2]


            return pdbs
        else:
            return -1

def elemToUnp(elem, autToUnp):
    arr = []
    for i in elem:
         arrElem = []
         for j in i:
             if j in autToUnp:
                if type(autToUnp[j]) == int:
                    arrElem.append(autToUnp[j])
                    if len(arrElem) == 2:
                        arr.append(arrElem)
    return arr

def getResList(resObj, pdbName, cat, eId):
    for chain in resObj[pdbName]['molecules']:

        if chain['entity_id'] == eId:
            for ch in chain['chains']:
               if ch['chain_id'] == cat:
                residues = {}
                for res in ch['residues']:
                    autRes = res['author_residue_number']
                    residues[res['residue_number']] = autRes
                return residues


def getResObj(pdbName):
    try:

            response = urllib.request.urlopen("https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/" + pdbName)
            data = json.load(response)
            return data
    except:
            print ('missing data from https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/'+ pdbName +'')
            return -1


def unpToAut(seq, residues, shift):

    unpAut = []
    i = 1
    while i <= len(seq):

        resNum = i - shift
        if resNum in residues:
           unpAut.append(residues[resNum])
           i += 1
        else:
           unpAut.append('-')
           i += 1
    return unpAut

def autToUnp(unpStart, unpEnd, residues, shift):
    autUnp = {}
    i = 1
    while i < len(residues):
        unpRes = i + shift
        if unpRes >= unpStart and unpRes <= unpEnd:
           autUnp[residues[i]] = unpRes
           i += 1
        else:
           autUnp[residues[i]] = 'u' + str(unpRes)
           i += 1
    return autUnp

def getPdbsEntities(pdb, ch):
    response = urllib.request.urlopen("http://repeatsdb.bio.unipd.it/ws/search?entry_type=repeat_region&id=" + pdb + ch + "&collection=repeat_region&show=ALL")
    data = json.load(response)
    if len(data) == 0:
        return [[],[],'']
    for elem in data:
        return [elem['units'], elem['insertions'], elem['classification']]




def unpsObjs(unps):

    fList = fInFolder()
    print(len(unipIds))
    print(len(fList))
    comp_list = [unp for unp in unps if unp not in fList]
    for unp in comp_list:
        print (unp)
        jsonObj = {}
        jsonObj['uniprot'] = []
        seq = getSeq(unp)
        if seq != -1:
           pdbs = getPdbs(unp, seq)
           if pdbs != -1:
               unipr = {}
               unipr['uniprotid'] = unp
               unipr['seq'] = seq
               unipr['pdbs'] = pdbs
               jsonObj['uniprot'].append(unipr)
               json_data = json.dumps(jsonObj)
               wrFile(json_data, unp)


def wrFile(jsondata, unp):
    f= open("files/" + unp + ".json","w+")
    f.write(jsondata)


def fInFolder():
    fList = os.listdir("files/")
    fTrimList = []
    for f in fList:
        f = f[:-5]
        fTrimList.append(f)
    return fTrimList



getId()
unpsObjs(unipIds)



