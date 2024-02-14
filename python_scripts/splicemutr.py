# The splicemutr python script library file

from binascii import hexlify, unhexlify
import numpy as np

#-----------------------------------------------------#

def calculate_similarity(counts_list):
    
    # calculate_similarity(): used in pep_similarity.py
    # Calculates the IEDB peptide similarity scores
    # inputs: 
    #  counts_list: list of alignment scores, one per kmer compared to the
    # returns:
    #  R: the peptide similarity score
    
    counts_np = np.array(counts_list)
    k = 4.87
    a = np.array(26)
    alpha = np.sum(np.exp(-1*k*(a - counts_list)))
    R = alpha/(1+alpha)
    return(float(R))

#-----------------------------------------------------#


class BytesIntEncoder:
    
    # The BytesIntEncoder class: used in extract_data.py
    # encodes the kmers as int and decodes the int to kmer

    @staticmethod
    def encode(b: bytes) -> int:
        return int(hexlify(b), 16) if b != b'' else 0

    @staticmethod
    def decode(i: int) -> int:
        return unhexlify('%x' % i) if i != 0 else b''
    
#-----------------------------------------------------#
    
def create_juncs(intron_dat):
    
    # create_juncs(): creates a string list of junctions
    # inputs:
    #  intron_dat: a pandas dataframe of intron_dat
    # outputs:
    #  juncs: the juncs string
    
    chr_dat = intron_dat.chr.tolist()
    start_dat = intron_dat.start.tolist()
    end_dat = intron_dat.end.tolist()
    juncs = [":".join([str(chr_dat[i]),str(start_dat[i]),str(end_dat[i])]) for i in range(len(chr_dat))]
    return(juncs)

#-----------------------------------------------------#

def create_specific_splicemutr(splice_dat, intron_file):
    # create_specific_splicemutr(): creating a splicemutr_dat file specific to one cancer type
    # inputs:
      # splice_dat: the full splicemutr output pandas dataframe
      # intron file: the leafcutter output for the current cancer type
    # outputs:
      # specific_splicemutr_file: the splice dat file that only contains rows found in the leafcutter output
    
    splice_dat_juncs = create_juncs(splice_dat)
    num_splice_dat = len(splice_dat_juncs)
    splice_dat["rows"] = range(num_splice_dat)
    intron_juncs = create_juncs(intron_file)
    splice_dat["juncs"] = splice_dat_juncs
    intron_file["juncs"] = intron_juncs
    
    intron_file = intron_file[intron_file.juncs.isin(splice_dat_juncs)]
    intron_juncs = intron_file.juncs.tolist()
    splice_dat = splice_dat[splice_dat.juncs.isin(intron_juncs)]
    splice_dat_juncs = splice_dat.juncs.tolist()
    intron_loc = [intron_juncs.index(j) for j in splice_dat_juncs]
    
    all_clusters = intron_file.clusterID.tolist()
    all_verdicts = intron_file.verdict.tolist()
    all_deltapsis = intron_file.deltapsi.tolist()

    splice_dat_new = pd.DataFrame({'cluster':[all_clusters[i] for i in intron_loc],
                      'chr': splice_dat.chr.tolist(),
                      'start': splice_dat.start.tolist(),
                      'end': splice_dat.end.tolist(),
                      'gene': splice_dat.gene.tolist(),
                      'tx_id': splice_dat.tx_id.tolist(),
                      'modified': splice_dat.modified.tolist(),
                      'is_UTR': splice_dat.is_UTR.tolist(),
                      'tx_junc_loc': splice_dat.tx_junc_loc.tolist(),
                      'pep_junc_loc': splice_dat.pep_junc_loc.tolist(),
                      'verdict':[all_verdicts[i] for i in intron_loc],
                      'deltapsi':[all_deltapsis[i] for i in intron_loc],
                      'error': splice_dat.error.tolist(),
                      'peptide': splice_dat.peptide.tolist(),
                      'tx_length': splice_dat.tx_length.tolist(), 
                      'start_stop': splice_dat.start_stop.tolist(),
                      'protein_coding': splice_dat.protein_coding.tolist(),
                      'rows': splice_dat.rows.tolist(),
                      'juncs': splice_dat.juncs.tolist()})
    
    return(splice_dat_new)

#-----------------------------------------------------#

def assign_kmers(genotypes_file,rows,hla_dir,hla_file):
    # assign_kmers(): assign kmers from mhc binding affinity data
    # inputs:
      # genotypes_file: the genotypes_file, including path
      # rows: the specific_splice_dat rows
      # hla_dir: the directory of the hla data
      # hla_file: the specific hla file without path
    # outputs:
      # kmers: imm. kmers for specific sample
      # sample_name: the specific sample name
        
    genotypes_file = genotypes_file.values.tolist()
    all_kmers = {}
    tumor_kmers = [set() for i in range(len(rows))]
    normal_kmers = [set() for i in range(len(rows))]
    kmers = [[] for i in range(len(rows))]
    #hlas = genotypes_file[0][0:6]
    #sample_type = genotypes_file[0][7]
    #sample_name = genotypes_file[0][8]
    hlas = genotypes_file[0][0].split(",")
    sample_name = genotypes_file[0][1]
    for hla in hlas:
        if (type(hla) is float):
            print("%s: break"%hla)
            continue
        else:
            with open(hla_file%(hla_dir,hla)) as geno_file:
                #print(hla_file%(hla_dir,hla))
                geno_list = geno_file.read().splitlines()
                geno_rows = [i.split('\t')[0] for i in geno_list]
                geno_kmers = [i.split('\t')[1] for i in geno_list]
                geno_dat = {int(geno_rows[i]):geno_kmers[i] for i in range(len(geno_rows))}
            for j in range(len(rows)):
                if rows[j] in geno_dat:
                    kmers[j].append(geno_dat[rows[j]])
        
    kmers = [":".join(k) for k in kmers]
    return(kmers,sample_name)

#-----------------------------------------------------#

def filter_kmers(kmers,groups_dframe):
    # inputs:
     # kmers: the kmers list
     # groups_dframe: the groups dataframe with columns groups, deltapsi, and rows
    # outputs:
     # kmers: the kmers with each row filtered by group and tumor or normal specificity
        
    groups_unique=list(set(groups_dframe.groups.tolist()))
    for group in groups_unique:
        groups_dframe_small=groups_dframe[groups_dframe["groups"]==group]
        normal_rows=groups_dframe_small[groups_dframe_small["deltapsi"]<0].rows.tolist() # this was originally >0
        tumor_rows=groups_dframe_small[groups_dframe_small["deltapsi"]>0].rows.tolist() # this was originally <0
        normal_kmers = set()
        [normal_kmers.update(kmers[i].split(":")) for i in normal_rows]
        tumor_kmers = set()
        [tumor_kmers.update(kmers[i].split(":")) for i in tumor_rows]
        for row in tumor_rows:
            kmers[row] = ":".join(list(set(kmers[row].split(":")).difference(normal_kmers)))
        for row in normal_rows:
            kmers[row] = ":".join(list(set(kmers[row].split(":")).difference(tumor_kmers)))
    return(kmers)

#-----------------------------------------------------#

def kmerize(pep,kmer_length):
    # kmerize() : kmerize a peptide
    # inputs:
     # pep: the peptide to be kmerized
     # kmer_length: the kmer length to kmerize the peptide over
    # outputs: 
     # the list of kmers from the peptide
    pep_len=len(pep)
    if len(pep) < kmer_length:
        pep += "Z"*(kmer_length-pep_len)
        pep_len=len(pep)
    return([pep[i:i+kmer_length] for i in range(pep_len-kmer_length+1)])
