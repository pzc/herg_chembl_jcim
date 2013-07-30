#! /usr/bin/env python

import urllib2
import json
import re

from sys import argv





codes = {'herg_human'  : 'Q12809' }



######

def looks_like_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False



def QueryChembl(accession=None):
    '''
    Query chembl
    '''
    ########################################################################

    # 1. Use UniProt accession to get target details

    print """
    # =========================================================
    # 1. Use UniProt accession to get target details
    # =========================================================
    """
    if not accession:
        accession = argv[1]
        #
    # test if we have a CHEMBL_TARGET_ID
    if accession.find("CHEMBL") == -1:     
      target_data = urllib2.urlopen("https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json" % accession).read()
      target_data =  json.loads(target_data)
      
      #print "Target Description: %s with %s" % (target_data['target']['description'],type(target_data['target']['description']))
      #print "Target CHEMBLID:    %s with %s" % (target_data['target']['chemblId'],type(target_data['target']['chemblId']))
      #
    #
    else:
      target_data = {}
      target_data['target'] = {}
      target_data['target']['chemblId'] = accession
    # 2. Get all bioactivties for target CHEMBL_ID

    print """

    # =========================================================
    # 2. Get all bioactivties for target CHEMBL_ID
    # =========================================================
    """

    bioactivity_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/targets/%s/bioactivities.json" % target_data['target']['chemblId']).read())

    print "Bioactivity Count:           %d" % len(bioactivity_data['bioactivities'])
    print "Bioactivity Count (IC50):    %d" % len([record for record in bioactivity_data['bioactivities'] if record['bioactivity_type'] == 'IC50'] )
    print "Bioactivity Count (Ki):      %d" % len([record for record in bioactivity_data['bioactivities'] if record['bioactivity_type'] == 'Ki' ])

   
    # 3. Get compounds with high binding affinity (IC50 < 100)
    
    print """

    # =========================================================
    # 3. Get compounds with binding affinity (IC50)
    # =========================================================
    """
    ic50_skip=0
    ki_skip=0
    inhb_skip=0

    count=0
    non_homo=0
    dr={}
    #for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('IC50', record['bioactivity_type']) and looks_like_number(record['value']) and float(record['value']) < 100]:
    #for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('IC50', record['bioactivity_type']) and looks_like_number(record['value']) ] :
    for bioactivity in [record for record in bioactivity_data['bioactivities'] if looks_like_number(record['value']) ] :
         
      #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
      #print
      if bioactivity['organism'] != 'Homo sapiens':
        non_homo+=1
        continue
      if re.search('IC50', bioactivity['bioactivity_type']):
        if bioactivity['units'] != 'nM':
          ic50_skip+=1
          continue
      elif re.search('Ki', bioactivity['bioactivity_type']):
        #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
        ki_skip+=1
        continue
      elif re.search('Inhibition', bioactivity['bioactivity_type']):
        inhb_skip+=1
      else:
        #print "Got activity type ", bioactivity['bioactivity_type']
        #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
        continue
        #
      #
      #print bioactivity['ingredient_cmpd_chemblid']
      cmpd_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/compounds/%s.json" % bioactivity['ingredient_cmpd_chemblid']).read())
      #print cmpd_data.items()
      my_smiles = cmpd_data['compound']['smiles']
      bioactivity['Smiles']=my_smiles
      dr[count] = bioactivity
      count+=1
      #
    #
    print "Skipped %i IC50 values" % ic50_skip
    print "Skipped %i Ki values" % ki_skip
    print "Skipped %i Inhibition values" % inhb_skip
    print "Skipped %i non-homo values" % non_homo
    #print dr
    import csv
    header = dr[0].keys()
    #print header
    
    ofile = "out_%s.csv" % (accession)
    dw = csv.DictWriter(open(ofile,'w'), delimiter='\t', fieldnames=header)
    dw.writerow(dict((fn,fn) for fn in header))
    for cpd in dr.keys():
        dw.writerow(dr[cpd])


    return
    # 4. Get assay details foe Ki actvity types


    print """

    # =========================================================
    # 4. Get assay details foe Ki actvity types
    # =========================================================
    """

    for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('Ki', record['bioactivity_type'], re.IGNORECASE)]:

      #print "Assay CHEMBLID: %s" % bioactivity['assay_chemblid']

      assay_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/assays/%s.json" % bioactivity['assay_chemblid']).read())

      #print "  %s" % assay_data['assay']['assayDescription']


if __name__ == '__main__':
    
    for name,accession in codes.items()[:]:
        print "Checking ", name
        QueryChembl(accession)
