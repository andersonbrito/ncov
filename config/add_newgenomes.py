#!/usr/bin/python
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with latest genomes from GISAID")
    parser.add_argument("--new-genomes", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--keep", required=True, help="TXT file with accession number of genomes to be included")
    parser.add_argument("--remove", required=True, help="TXT file with accession number of genomes to be removed")
    parser.add_argument("--output", required=True, help="FASTA file containing filtered sequences")
    args = parser.parse_args()

    genomes = args.genomes
    new_genomes = args.new_genomes
    keep = args.keep
    remove = args.remove
    outfile = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/run7_20210204_variants/ncov/'
    # genomes = path + "pre-analyses/gisaid_hcov-19.fasta"
    # new_genomes = path + "pre-analyses/new_genomes.fasta"
    # keep = path + 'config/keep.txt'
    # remove = path + "config/remove.txt"
    # outfile = path + "pre-analyses/temp_sequences.fasta"


    # store only new sequences in a dictionary, ignoring existing ones
    print('\n### Loading new sequences\n')
    newly_sequenced = {}
    for fasta in SeqIO.parse(open(new_genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        print(id)
        if id not in newly_sequenced.keys(): # avoid potential duplicates
            newly_sequenced[id] = str(seq)
    # print(newly_sequenced)

    # create a list of the existing sequences
    print('\n### Loading GISAID dataset\n')
    # all_sequences = {}
    all_sequences = []
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        # id, seq = fasta.description, fasta.seq
        id = fasta.description
        # print(id)
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        # if id not in newly_sequenced.keys(): # avoid potential duplicates
        if id not in newly_sequenced:  # avoid potential duplicates
            if "Yale-" not in id: # change this line to match you lab's unique identifier
                all_sequences.append(id)
                # all_sequences[id] = str(seq)

    # create a list of sequences to be added in all instances
    # keep_sequences = {}
    keep_sequences = []
    mismatch = []
    for id in sorted(open(keep, "r").readlines()):
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id not in newly_sequenced:
                if id not in keep_sequences:
                    if 'Yale-' not in id: # change this line to match you lab's unique identifier
                        try:
                            # keep_sequences[id] = all_sequences[id]
                            keep_sequences.append(id)
                        except:
                            mismatch.append(id)

    # create a list of sequences to be ignored in all instances
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)


    # export only sequences to be used in the nextstrain build
    c = 1
    # sequences = {**keep_sequences, **newly_sequenced}
    # sequences = keep_sequences
    print('\n### Exporting sequences\n')
    exported = []
    with open(outfile, 'w') as output:
        for fasta in SeqIO.parse(open(genomes), 'fasta'):
            id, seq = fasta.description, fasta.seq
            if id not in remove_sequences:
                if id in keep_sequences: # filter out unwanted sequences
                    entry = ">" + id + "\n" + str(seq).upper() + "\n"
                    exported.append(id)
                    output.write(entry)
                    # if 'Yale-' in id: # change this line to match you lab's unique identifier
                    #     print('* ' + str(c) + '. ' + id)
                    # else:
                    print(str(c) + '. ' + id)
                    c += 1
            # else:
            #     c -= 1
        for id, seq in newly_sequenced.items():
            print('* ' + str(c) + '. ' + id)
            entry = ">" + id + "\n" + seq.upper() + "\n"
            exported.append(id)
            output.write(entry)
            c += 1



    # mismatched sequence headers
    if len(mismatch) > 0:
        print('\n### Possible sequence header mismatches\n')
        m = 1
        for id in mismatch:
            print(str(m) + '. ' + id)
            m += 1
    else:
        print('\nNo sequence name mismatches found...')


    # # excluding sequences
    # print('\n### Excluding sequences ###\n')
    # e = 1
    # for id in remove_sequences:
    #     print(str(e) + '. ' + id)
    #     e += 1

    print('\n### Final result\n')

    print('Lab file contains ' + str(len(newly_sequenced)) + ' sequences')
    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(newly_sequenced)) + ' newly sequenced genomes were added')
    print(str(len(remove_sequences)) + ' genomes were REMOVED according to remove.txt\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')