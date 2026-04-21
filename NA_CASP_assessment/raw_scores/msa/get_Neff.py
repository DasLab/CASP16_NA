import numpy as np
from glob import glob
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse

def get_Neffs(msa_file, thresholds,msa_format='fasta',method='af2',per_residue=False):
    # read in the msa
    sequences = [record.seq for record in SeqIO.parse(msa_file,msa_format)]
    num_sequences =len(sequences)

    if per_residue:
        per_residue_weights = {}
        for threshold in thresholds:
            per_residue_weights[threshold] = np.array([np.zeros(len(sequences[0])) for x in range(num_sequences)])

    gaps = ['-','.']

    # Initialize weights
    weights = {}
    for threshold in thresholds:
        weights[threshold] = np.zeros(num_sequences)
    #print('MSA contains: ',num_sequences)
    
    # Calculate weights
    for i in tqdm(range(num_sequences)):
        # just shorten the run time a little
        for j in range(i,num_sequences):
            if method == 'af2':
                # Ref 76 from the AF2 paper explicitly says to keep the self score
                
                # RCK for sequence identity we only pay attention to regions where at least one is not gap
                # it seems from AF2 methods that this is pairwsie so should be done in this loop instead of the previous remove insert functionality
                nongapped_index = [n for n in range(len(sequences[i])) if sequences[i][n] not in gaps or sequences[j][n] not in gaps]
                #gapped_index = [n for n in range(len(sequences[i])) if sequences[i][n] in gaps and sequences[j][n] in gaps]                
                seqA = "".join([sequences[i][index] for index in nongapped_index])
                seqB = "".join([sequences[j][index] for index in nongapped_index])
                #seqA_gapped = "".join([sequences[i][index] for index in gapped_index])
                #seqB_gapped = "".join([sequences[j][index] for index in gapped_index])
                problem_chars = [x for x in seqA+seqB if x not in ['A','C','G','U']+gaps]
                if len(problem_chars)>0:
                    print('problem',problem_chars)
                    
                # Calculate sequence identity between sequences i and j
                identity = sum(res1 == res2 for res1, res2 in zip(seqA,seqB)) / len(nongapped_index)
                # identity is normalized identically
                
                for threshold in thresholds:
                    if identity > threshold:
                        weights[threshold][i] += 1
                        # for shorten run time only calculate each pair once
                        weights[threshold][j] += 1

                        if per_residue:
                            non_gap_A = [n for n in range(len(sequences[i])) if sequences[i][n] not in gaps]
                            for index in non_gap_A:
                                per_residue_weights[threshold][i][index] += 1
                                
                            non_gap_B = [n for n in range(len(sequences[j])) if sequences[j][n] not in gaps]
                            for index in non_gap_B:
                                per_residue_weights[threshold][j][index] += 1
                    
    
            elif method == 'neffy_symmetric':
                seqA = sequences[i]
                seqB = sequences[j]
                    
                # Calculate sequence identity between sequences i and j
                identity = sum(res1 == res2 for res1, res2 in zip(seqA,seqB)) / len(seqA)
    
                for threshold in thresholds:
                    if identity > threshold:
                        weights[threshold][i] += 1
                        # for shorten run time only calculate each pair once
                        weights[threshold][j] += 1
            elif method == 'neffy_non_symmetric':                
                # RCK for sequence identity we only pay attention to regions where at least one is not gap
                # it seems from AF2 methods that this is pairwsie so should be done in this loop instead of the previous remove insert functionality
                nongapped_index = [n for n in range(len(sequences[i])) if sequences[i][n] not in gaps or sequences[j][n] not in gaps]
                #gapped_index = [n for n in range(len(sequences[i])) if sequences[i][n] in gaps and sequences[j][n] in gaps]
                
                seqA = "".join([sequences[i][index] for index in nongapped_index])
                seqB = "".join([sequences[j][index] for index in nongapped_index])
                #seqA_gapped = "".join([sequences[i][index] for index in gapped_index])
                #seqB_gapped = "".join([sequences[j][index] for index in gapped_index])
                problem_chars = [x for x in seqA+seqB if x not in ['A','C','G','U']+gaps]
                if len(problem_chars)>0:
                    print('problem',problem_chars)
                    
                # the emulate Neffy
                total = sum(res1 == res2 for res1, res2 in zip(seqA,seqB))
                non_gap_A = len([n for n in range(len(sequences[i])) if sequences[i][n] not in gaps])
                non_gap_B = len([n for n in range(len(sequences[j])) if sequences[j][n] not in gaps])
                # identity is normalized deifferently for each
    
                for threshold in thresholds:
                    # the emulate Neffy # identity is normalized 
                    if total/non_gap_A > threshold:
                        weights[threshold][i] += 1
                    if total/non_gap_B > threshold:
                        weights[threshold][j] += 1

            else:
                print('ERROR: only af2, neffy_non_symmetric, neffy_symmetric implemented')

        # per residue sums the weight that have a residue at that position
        # eg I get NxN matrix with weights from the above
        # for Neff it is taking each row
        # for per residue I only add for nongapped for seqA or seqB?
        
    # Calculate the weight for sequence i
    for threshold in thresholds:
        # i==j get counted twice, so -1 to correct
        weights[threshold] = 1 / (weights[threshold]-1)
        if per_residue:
            non_zero = per_residue_weights[threshold] != 0
            per_residue_weights[threshold][non_zero] = 1 / (per_residue_weights[threshold][non_zero]-1)

    # Calculate Neff
    Neff = {}
    per_residue_Neff = {}
    for threshold in thresholds:
        Neff[threshold] = np.sum(weights[threshold])
        #print(f'Neff at {threshold} = {Neff[threshold]}')
        if per_residue:
            per_residue_Neff[threshold] = per_residue_weights[threshold].sum(axis=0)
            #print(f'Per residue Neff at {threshold} = {per_residue_Neff[threshold]}')
        #print(weights[threshold])
    if per_residue:
        return num_sequences, Neff, per_residue_Neff
    else:
        return num_sequences, Neff

# python get_Neff.py --msa_file_pattern "from_kaggle/R*fasta" --thresholds 0.62 0.8 --per_residue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate effective number of sequences (Neffs).")
    parser.add_argument("--msa_file_pattern", help="Path to the msa files")
    parser.add_argument("--thresholds", nargs="+", type=float, help="List of thresholds (space-separated floats)")
    parser.add_argument("--msa_format", default="fasta", help="Format of the MSA file (default: fasta)")
    parser.add_argument("--method", default="af2", help="Method used for Neff calculation (default: af2)")
    parser.add_argument('--per_residue', action='store_true')
    parser.add_argument("--output_prefix", default="neff_results", help="Prefix for output CSV files")

    args = parser.parse_args()

    msa_files = glob(args.msa_file_pattern)
    all_neffs = {}  # Dictionary to store Neffs for all files
    all_per_residue_neffs = {} # Dictionary to store the per residue for all files

    for msa_file in msa_files:
        print(f"Processing file: {msa_file}")
        try:
            results = get_Neffs(msa_file, args.thresholds, args.msa_format, args.method,args.per_residue)
            if args.per_residue:
                N, neffs, per_residue_neffs = results
                all_per_residue_neffs[msa_file] = per_residue_neffs

            else:
                N, neffs = results
            neffs['N'] = N
            all_neffs[msa_file] = neffs

        except Exception as e:
            print(f"Error processing {msa_file}: {e}")

        # Save Neffs to CSV
        neff_data = {}
        for threshold in args.thresholds + ['N']:
            neff_data[threshold] = []
        for msa_file, neff_results in all_neffs.items():
            for threshold, neff_value in neff_results.items():
                neff_data[threshold].append(neff_value)
        # make the dataframe
        df_neffs = pd.DataFrame(neff_data,index=all_neffs.keys())
        neff_output_file = f"{args.output_prefix}_neffs.csv"
        df_neffs.to_csv(neff_output_file)
        print(f"Neff results saved to: {neff_output_file}")

        # Save per-residue Neffs to CSV
        for threshold in args.thresholds:
            all_residue_data = {}
            for msa_file, per_residue_neff_results in all_per_residue_neffs.items():
                all_residue_data[msa_file] = per_residue_neff_results[threshold]
            df_per_residue = pd.DataFrame.from_dict(all_residue_data,orient='index') #Orient index is a cool trick
            per_residue_output_file = f"{args.output_prefix}_per_residue_neffs_{threshold}.csv"
            df_per_residue.to_csv(per_residue_output_file)
            print(f"Per-residue Neff results (threshold={threshold}) saved to: {per_residue_output_file}")




'''
The MSA depth analysis was based on computing the normalized number of effective sequences (Neff) for each position of a query sequence. 
Per-residue Neff values were obtained by counting the number of non-gap residues in the MSA for this position and weighting the sequences 
using the Neff scheme76 with a threshold of 80% sequence identity measured on the region that is non-gap in either sequence.

'''
# do we double count??? j>i? --- yes I check ref 76 from AF2 and they did indeed count for all j including i=j!

# are we meant to drop duplicates, does that change to scoring? -- no it does not appear they drop this, in fact the explictly call out an algorithmn that drops at X% seq id as a difference
# to be fair, unlike the Neff reference, AF2 paper did seem to pull MSA from databases already void of redudancy so in principle they were doing this, but not in the Neff calculation
# either way score change is minimal

# essentially we are using the Neff from 76 but seqid from Af2 -- should just calc both 80 (AF2) and 62 (original Neff cutoff)

# Neffy differs in the the identity calculation defines identity per seqA or B based on their number of nongap.... 
