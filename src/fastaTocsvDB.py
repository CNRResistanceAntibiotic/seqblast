from Bio import SeqIO
import os

def fastaRead(filename):
    recList = []
    for rec in SeqIO.parse(open(filename, 'r'),'fasta'):
        recList.append(rec)
    print '\nRead %i record(s)\n' % len(recList)
    return recList

def recParse(recList):
    prbDict = {}
    recDict = {}
    fctDict={'Bla':'Beta-lactam','AGly':'Aminoglycoside','Iso':'Isoniazid','Flq':'Quinolone','Plx':'Colistin',
             'Tet':'Cycline','Efl':'Efflux','Fcyn':'Fosfomycin','Eta':'Ethambutol','Reg':'Efflux','Efa':'Elfamycin',
             'Eti':'Etionamide','Fus':'Fusidic_acid','Gly':'Glycopeptide','Nit':'Nitro-imidazole',
             'Lpp':'Lipopeptide','Ozd':'Oxazolidinone|Phenicol','Phe':'Phenicol','Pyr':'Pyrazinamide',
             'Rif':'Rifampin','MLS':'Lincosamide|Macrolide|Streptogramin','Sul':'Sulfonamide','Tmt':'Trimethoprim',
             'AGlyFlqn':'Aminoglycoside|Quinolone','Tun':'Divers','Bac':'Divers', 'Ble':'Divers', 'Imp':'Permeability',
             'ATB':'Divers'}

    for rec in recList:
        grp_nb, grp_name, id_name, id_nb = rec.id.split('__')
        try :
            entry_name = rec.description.split(';')[2]
            function_grp_names = rec.description.split(';')[3]
            if function_grp_names == '':
                function_grp_names = grp_name.split('_')[1]
            function_grp_names = fctDict[function_grp_names]
            dnaAcc = rec.description.split(';')[4]
        except IndexError:
            entry_name = id_name
            function_grp_names = grp_name.split('_')[1]
            function_grp_names = fctDict[function_grp_names]
            dnaAcc = rec.description.split('_')[-1]

        if entry_name != id_name:
            entry_name = id_name

        dna_sequence = rec.seq
        try :
            prot_sequence = rec.seq.translate(table="Bacterial", cds=True, to_stop=True)
        except:
            prot_sequence = ''

        if prot_sequence == '' or '*' in prot_sequence:
            prbDict[rec.id] = {'entry_name':entry_name,'dna_sequence':dna_sequence,'dna_accession':dnaAcc,
                           'prot_sequence':prot_sequence,'function_grp_names':function_grp_names}
        else:
            recDict[rec.id] = {'entry_name':entry_name,'dna_sequence':dna_sequence,'dna_accession':dnaAcc,
                           'prot_sequence':prot_sequence,'function_grp_names':function_grp_names}
        print '\n%s: entry_name:%s function_grp_names:%s dna_accession:%s' % (rec.id, entry_name, function_grp_names, dnaAcc)
        print 'dna : %s' % dna_sequence
        print 'prot: %s' % prot_sequence
        print ''
    return recDict, prbDict

def csvWrite(recList, Dict, filename):
    header = ['entry_name', 'function_grp_names', 'dna_accession', 'prot_sequence', 'dna_sequence']
    f = open(filename, 'w')
    f.write('key\t'+'\t'.join(header)+'\n')
    for rec in recList:
        line = rec.id
        if rec.id in Dict.keys():
            for item in header:
                line = line + '\t' + Dict[rec.id][item]
            line = str(line + '\n')
            f.write(line)
    f.close()

fnaFile = 'cnr/db_arm_5.fna'
outfile = os.path.splitext(fnaFile)[0]+'.csv'
prbfile = os.path.splitext(fnaFile)[0]+'_prb.csv'
recList = fastaRead(fnaFile)
recDict, prbDict = recParse(recList)
csvWrite(recList, recDict, outfile)
csvWrite(recList, prbDict, prbfile)
