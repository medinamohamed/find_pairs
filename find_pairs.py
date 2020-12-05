import sys
import pandas as pd 
import time

import numpy as np

import sys
from sys import stdout


def main():

    #1. Read the data from the PAF file
    paf = read_paf_file(sys.argv[1])

    #2. Read the data from the summary file 
    summ = read_summ_file(sys.argv[2])

    #3. List with all the reads in common in both files
    reads_list = merge_lists(summ,paf) # param

    #4. Create a new dataframe with needed attributes
    df = merge_dataframe(reads_list,summ,paf)

    #5. List of all the reads that are in a pair for a given chromosome
    reads = reads_list_per_chrm(df,'chr1') #test for chrom. 12

    #6. Filter the dataframe with all the reads in a pair.
    df = df[df.read_id.isin(reads)] 

    #7.Save file in tsv format
    df.to_csv(sys.argv[3], sep='\t')

def read_paf_file(paf_file):
   
    paf = pd.read_table(sys.argv[1]) 

    paf.drop(paf.iloc[:,12:18],inplace=True,axis=1) 
    paf.columns = ["read_id","query_length", "query_start","query_end","sign", "chr", "target_length","target_start","target_end",
                "bases","bases_gaps","mapq"]

    # Generates a sorted list of the reads in the PAF file
    paf = paf.sort_values(by ='read_id')
    paf = paf.reset_index(drop=True)

    # Remove reads that have a quality lower than 60 
    paf = paf[paf.mapq >= 60]

    return paf

def read_summ_file(summ_file):

    summ = pd.read_table(sys.argv[2])

    # Generate a sorted list of the reads in the summary file
    summ = summ.sort_values(by ='read_id')
    summ = summ.reset_index(drop=True)

    return summ


def merge_lists(summ,paf):
    # Returns a list of all the reads that are both in the paf file and summary file
    summ_ids = summ['read_id']
    paf_ids  = paf['read_id']

    return np.intersect1d(summ_ids,paf_ids)


def merge_dataframe(reads_list,summ,paf):

    # Reduce PAF AND SUMM file to the ids in commons
    paf = paf[paf.read_id.isin(reads_list)]
    summ = summ[summ.read_id.isin(reads_list)]

    # Keep needed  attributes in PAF and summ
    paf = paf.reindex(['read_id','target_start','target_end','target_length','sign','chr','mapq'],axis=1)

    summ = summ.reindex(['read_id','channel','start_time','duration'],axis=1)

    # with ids in common merge dataframe 
    df = summ.set_index('read_id').join(paf.set_index('read_id')).reset_index()

    return df
    
def filter_reads(df, chrm, pore_id):
    # returns a dataframe for one pore and one chromosome
    df = df.loc[(df.channel == pore_id) & (df.chr == chrm)]
  
    return df

def filter_strands(df,sign):

    if (sign == '+'):
        df = df.loc[df.sign == '+']

    if (sign == '-'):
        df = df.loc[df.sign == '-']

    return df   

def find_overlap(df,first_read_id, sec_read_id,start_first_read, end_first_read,sign_first_read):
    
    start_sec_read = df[df['read_id']== sec_read_id]['target_start'].values[0]
    end_sec_read = df[df['read_id']== sec_read_id]['target_end'].values[0]
    
    x = range(start_first_read,end_first_read)
    y = range(start_sec_read,  end_sec_read)
    
    sign_sec_read = df[df['read_id']== sec_read_id]['sign'].values[0]
    
    
    if  (len(set(x) & set(y)) != 0) and (sign_first_read != sign_sec_read  ):
        return True
    else :
        return False

def find_one_pair(read_id,df,window): 

    start_first_read = df[df.read_id == read_id]['target_start'].values[0]
    end_first_read  = df[df.read_id == read_id]['target_end'].values[0]

    sign = df[df['read_id']== read_id]['sign'].values[0]

    # df at this point is already filtered by chrm & pore
    df = filter_strands(df,'-')

  
    df = df.loc[(start_first_read  < df.target_start + window) & (start_first_read  > df.target_start  -window)]           
    reads_ids = df['read_id']

    
    for i in reads_ids:

        if (find_overlap(df,read_id,str(i),start_first_read,end_first_read,sign)): 
            # find overlap version end 
            return [read_id, str(i)]       

def find_other_pairs(df): 

    #df at this point is already filtered by chrm & pore
    df_pos = filter_strands(df,'+')

    reads_ids = df_pos['read_id']

    all_pairs = []
    for read_id in reads_ids:
        if(find_one_pair(read_id,df,150)):
            all_pairs.append(find_one_pair(read_id,df,150))
            
    return all_pairs    

def reads_list_per_chrm(df,chrm):
 
    pores = np.array((df.channel.unique()))

    reads_list_per_chrm = []

    for i in range(0,len(pores)):

        stdout.write("\r%d" % ((i/len(pores))*101)+ " %\t")
        stdout.flush()
        
        df_filtered = filter_reads(df,chrm,pores[i])
        tab = find_other_pairs(df_filtered)
        

        for a,b in tab:
            reads_list_per_chrm.append(a)
            reads_list_per_chrm.append(b)

  
    return reads_list_per_chrm

if __name__ == "__main__":
    main()
