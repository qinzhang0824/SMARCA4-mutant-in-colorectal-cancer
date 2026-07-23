import os
import re
import sys
import click
import argparse
import math
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(description="TMB",formatter_class=RawTextHelpFormatter)
parser.add_argument('--mutinfo',help="high",required=True)
#parser.add_argument('--fracinfo',help="high",required=True)
#parser.add_argument('--cell_info',help="",required=True)
#parser.add_argument('--high_cutoff',help="high",required=True)
parser.add_argument('--fpkm',help="",required=True)
parser.add_argument('--dir',help="low",required=True)
argv = parser.parse_args()



cell_index={}

cell_sample={}

sample_mut={}
sample_wild={}
sample_mut_status={}

with open(argv.mutinfo,'r') as ff:
    for each in ff:
        if each.startswith('Sample'):
            head=each.rstrip().split('\t')
            #stage_index=head.index("Stage")
            pat_index=head.index("Sample")
            group_index=head.index("SMARCA4")
            continue
        else:
            info=each.rstrip().split('\t')
            pat_info=info[pat_index].split('-')
            #pat='-'.join(pat_info[0:4])
            pat=info[pat_index]
            sample_mut_status[pat]=info[group_index]
            #if info[stage_index]=="III" or info[stage_index]=="IV":
            if info[group_index]=="Mut":
                sample_mut[pat]=info[group_index]
            else:
                sample_wild[pat]=info[group_index]
print(len(sample_mut),len(sample_wild))
sample_index=[]
sample_order=[]
target_sample={}

fpkm_outfile=open(os.path.join(argv.dir,'tcga_normalized_counts_for_GSEA_input.gct'),'w')
with open(argv.fpkm,'r') as ff:
    for each in ff:
        tmp=[]
        if each.startswith('#'):
            fpkm_outfile.write(each)
        elif each.startswith('gene_id'):
            head=each.rstrip().split('\t')
            #fpkm_outfile.write(head[0]+'\t'+"gene_name"+'\t'+'\t'.join(target_sample)+"\n")
            for i in head[1:]:
                name=i.split('-')
                samid=re.search('\d*',name[3])
                patient_post=samid.group()
                tar='-'.join(name[0:3])+"-"+patient_post
                if tar in sample_mut:
                    target_sample[head.index(i)]=tar
                    sample_index.append(head.index(i))
                    sample_order.append(tar)
            for j in head[1:]:
                nam=j.split('-')
                pp=re.search('\d*',nam[3])
                pp_pos=pp.group()
                ta='-'.join(nam[0:3])+"-"+pp_pos
                if ta in sample_wild:
                    target_sample[head.index(j)]=ta
                    sample_index.append(head.index(j))
                    sample_order.append(ta)
            #print len(sample_index)
            fpkm_outfile.write("NAME"+'\t'+"Description"+'\t'+'\t'.join(sample_order)+"\n")
        else:
            info=each.rstrip().split('\t')
            for each_index in sample_index:
                if info[each_index]!="NA":
                    tmp.append(info[each_index])
            if info[0].split('|')[0]=="?" or len(tmp) != len(sample_index):
                continue
            else:
                fpkm_outfile.write(info[0].split('|')[0]+'\t'+"naa"+'\t'+'\t'.join(tmp)+"\n") 


out=os.path.join(argv.dir,'tcga_sample_group_for_GSEA_input.cls')
outfile=open(out,'w')

outfile.write(str(len(sample_order))+"\t"+"2"+"\t"+"1"+"\n")
outfile.write("#"+"\t"+"Mut"+"\t"+"Wildtype"+"\n")

mut_status=[]
for each_sam in sample_order:
    mut_status.append(sample_mut_status[each_sam])
outfile.write("\t".join(mut_status)+"\n")
outfile.close()

