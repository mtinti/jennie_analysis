
#separate by random seq ''.join([random.choice('AGTC') for x in range(10)]) 

from Bio import SeqIO
import random
from string import strip
import pandas as pd


def make_dict_from_fasta(in_fasta):
    res = {}
    for record in SeqIO.parse("vsgs_all_cds.fasta", "fasta"):
        #print '|'.join(record.description.split('|')[3:])
        #print record.id
        if '|' in record.description:
            desc = '|'.join(record.description.split('|')[3:])
        else:
            desc=record.description
        
        res[record.id]=strip(desc)
    return res


def make_dict_from_gff(in_gff):
    res = {}
    def make_dict(in_string):
        res = {}
        for item in in_string.replace('%2C',' ').split(';'):
            res[strip(item.split('=')[0])]= strip(item.split('=')[1])
        return res
    for line in open(in_gff):
        if ('ID=' in line) and ('Parent' not in line):
            item_list=line.split('\t')
            temp = make_dict(item_list[8])
            #print temp
            res[temp['ID']]=temp['description']
    return res
            
        

if __name__ == '__main__':
    in_gff = 'TriTrypDB-33_TbruceiTREU927.gff'
    dict_1 = make_dict_from_gff(in_gff)
    in_fasta = 'assemble_vsg_genome/vsgs_all_cds.fasta'
    dict_2 = make_dict_from_fasta(in_fasta)
    dict_1.update(dict_2)
    print dict_1['Tb427VSG-7565']
    df = pd.DataFrame.from_csv('ordered_results_deseq_2.csv')
    df['desc']=[dict_1[n] for n in df.index.values]
    print df.head()
    df.to_csv('final_table.csv')

'''
def insert_newlines(string, every=60):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

    
#concatenate the vsg sequences in one long fasta seq    
handle = open('concatenated_vsg.fasta','w')
header = '>fake_vsg_chr | organism=Trypanosoma_brucei | version=May 2017 | length={n} | SO=chromosome\n'


start=1
end=1
fasta_seq = ''.join([random.choice('AGTC') for x in range(20)]) 
for record in SeqIO.parse("vsgs_all_cds.fasta", "fasta"):
    if 'Tb427' in record.id:
        fasta_seq+=str(record.seq)
        fasta_seq+=''.join([random.choice('AGTC') for x in range(20)])
header=header.replace('{n}',str(len(fasta_seq)))
handle.write(header+insert_newlines(fasta_seq, every=60))
handle.close()


#add gtf information for the vsg, only if not present in other chr
all_seq = ''
for record in SeqIO.parse("TriTrypDB-33_TbruceiTREU927_Genome.fasta", "fasta"):
    all_seq+=str(record.seq)

handle = open('add.gtf','w')
last_field = 'transcript_id "{gene_id}:mRNA"; gene_id "{gene_id}"';
template_string = 'fake_vsg_chr'+'\t'+'tryps.rockefeller.edu'+'\t'+'{feature}'+'\t'+'{start}'+'\t'+'{end}'+'\t'+'.'+'\t'+'+'+'\t'+'{frame}'+'\t'+last_field

assembled = ''
a=0
b=0
for record in SeqIO.parse("vsgs_all_cds.fasta", "fasta"):
    a+=1
    #print a
    temp_seq = str(record.seq)
    print a,record.id
    if 'Tb427' in record.id:
        if temp_seq not in assembled:
            if temp_seq not in all_seq:
                print 'found' 
                start = fasta_seq.find(temp_seq)
                end = start+len(temp_seq)
                temp = template_string.replace('{feature}','exon').replace('{start}',str(start+1)).replace('{end}',str(end)).replace('{frame}','.').replace('{gene_id}',record.id)
                handle.write(temp+'\n')
                temp = template_string.replace('{feature}','CDS').replace('{start}',str(start+1)).replace('{end}',str(end)).replace('{frame}','.').replace('{gene_id}',record.id)
                handle.write(temp+'\n')
                b+=1
                assembled+=temp_seq
print a,'parsed'
print b,'added'
handle.close()
'''    





