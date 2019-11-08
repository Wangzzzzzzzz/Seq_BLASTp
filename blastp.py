from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import os
import numpy as np 
import pandas as pd
import json

# extract the sequences based on the dataset we have 
def obtain_seq(score, sequences):
    score  = './dataFile/' + score 
    sequences = './dataFile/' + sequences
    score_file = pd.read_csv(score, sep='\t')
    wt_1 = score_file.iloc[:,0]
    wt_2 = score_file.iloc[:,1]
    # concatinate the wt_1 and wt_2 file 
    wild_type = list(wt_1) + list(wt_2)
    # get only the unique items and then sort them
    wild_type = set(wild_type)
    
    # read in the sequence file and compile into a dictionary
    seq_dict = {}
    with open(sequences, 'r') as s:
        for l in s:
            # consider only the wildtype
            content = l.split('\t')
            sequence = content[1].replace(' ','').replace('\n','')
            header = content[0]
            header_info = header.split('_')
            if (len(header_info) == 2) and (header in wild_type):
                seq_dict[header] = sequence

    return seq_dict

def write_to_fasta(seq_dict, fasta_name):
    # generate the SeqRecord 
    seqRec_list = [SeqRecord(Seq(seq_str, IUPAC.protein), id=seq_key)
                   for seq_key, seq_str in seq_dict.items()]
    with open(fasta_name, 'w') as output_hl:
        SeqIO.write(seqRec_list, output_hl, "fasta")

def run_blast(query, subject):
    output = NcbiblastpCommandline(
        query=query, subject=subject, outfmt=5, evalue=0.05)()[0]
    # read all parser result
    blast_result_record = list(NCBIXML.parse(StringIO(output)))
    return blast_result_record

def output_json(blast_res, output_name):
    #output the result in JSON file, compile the results into a conplex dictionary
    dict_blast_res = {}
    # each object in the list is the alignments find for one
    # one sequence in the query set
    for query_seq in blast_res:
        # query seq's name
        query_seq_name = query_seq.query[0:6]
        if len(query_seq.alignments):
            # create an empty dictionary for that query seq & add this query seq to dictionary
            dict_blast_res[query_seq_name] = dict()
        else:
            continue
        # go through all alignments of all one sequence in query
        for alignment in query_seq.alignments:
            # count of hsp, since it does not have unique identifier
            hsp_id = 1
            # subject seq (aligned sequence in database)'s name
            subjt_seq_name = alignment.title[0:6]
            if len(alignment.hsps):
                # add the subjt name to the sub dict, and create an empty dict for hsp
                dict_blast_res[query_seq_name][subjt_seq_name] = dict()
            else:
                continue
            # go through all of the high score pairs (represent the regions in the hit seq
            # that contains significant alignment, usually there is only one)
            for hsp in alignment.hsps:
                hsp_discriptor = {
                    "score": hsp.score,
                    "E_Value": hsp.expect,
                    "Perc_Identity": hsp.identities/alignment.length*100,
                    "Num_of_Identity": hsp.identities,
                    "Aligned_length": hsp.align_length,
                    "Align_Seq_Length": alignment.length,
                    "Query_Seq_Length": query_seq.query_length,
                    "query_start_pos": hsp.query_start,
                    "aligned_query_seq": hsp.query,
                    "aligned_match_pts": hsp.match,
                    "aligned_sbjct_seq": hsp.sbjct,
                    "sbjct_start_pos":hsp.sbjct_start,
                }
                # filter of the Perc_Identity, change if needed
                if hsp_discriptor['Perc_Identity'] < 25:
                    continue
                # add the hsp_id to dict & attach the discriptor of hsp to the dict
                dict_blast_res[query_seq_name][subjt_seq_name][hsp_id] = hsp_discriptor
                hsp_id += 1

            # remove empty dicts
            if not dict_blast_res[query_seq_name][subjt_seq_name]:
                del dict_blast_res[query_seq_name][subjt_seq_name]

        if not dict_blast_res[query_seq_name]:
            del dict_blast_res[query_seq_name]

    #if there are something to output
    if dict_blast_res:
        with open('./blast_rs/' + output_name, 'w') as output_hl:
            json.dump(dict_blast_res, output_hl,indent=4)
        return 1
    
    # if there is nothing to output:
    else:
        return 0


def main():
    if not os.path.exists('./fasta_db'):
        os.mkdir('./fasta_db')

    if not os.path.exists('./blast_rs'):
        os.mkdir('./blast_rs')

    # run blast on skempi_v1 and skempi_v2
    skempi_v1 = obtain_seq('./SKP1402m.ddg.txt','./SKP1402m.seq.txt')
    write_to_fasta(skempi_v1, './fasta_db/skempi_v1.fasta')
    skempi_v2 = obtain_seq('./3G_S487_test_dataset_top1.txt','./skempi_v2.singlemut.mut4.seq.txt')
    write_to_fasta(skempi_v2, './fasta_db/skempi_v2.fasta')
    blast_skempi_v1_skempi_v2 = run_blast(query='./fasta_db/skempi_v1.fasta',
                                          subject='./fasta_db/skempi_v2.fasta')
    blast_skempi_v1_skempi_v1 = run_blast(query='./fasta_db/skempi_v1.fasta',
                                          subject='./fasta_db/skempi_v1.fasta')
    blast_skempi_v2_skempi_v2 = run_blast(query='./fasta_db/skempi_v2.fasta',
                                          subject='./fasta_db/skempi_v2.fasta')

    # run blast on skempi_v1 and NM
    NM = obtain_seq('./NM_test.scores.txt','NM_test.seq.txt')
    write_to_fasta(NM, './fasta_db/NM.fasta')
    blast_skempi_v1_NM = run_blast(query='./fasta_db/skempi_v1.fasta',subject='./fasta_db/NM.fasta')


    # output the json file for all of the above results
    if not output_json(blast_skempi_v1_skempi_v1, 'skempi_v1_skempi_v1.json'):
        print('No alignments for skempi_v1 & skempi_v1')
    if not output_json(blast_skempi_v1_skempi_v2, 'skempi_v1_skempi_v2.json'):
        print('No alignments for skempi_v1 & skempi_v2')
    if not output_json(blast_skempi_v2_skempi_v2, 'skempi_v2_skempi_v2.json'):
        print('No alignments for skempi_v2 & skempi_v2')
    if not output_json(blast_skempi_v1_NM, 'skempi_v1_NM.json'):
        print('No alignments for skempi_v1 & NM')

if __name__ == "__main__":
    main()

