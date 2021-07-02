#!/usr/bin/python
import argparse
from Bio import SeqIO
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA or JSON file with latest genomes from GISAID")
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

    # # genomes = path + "pre-analyses/provision.json"
    # genomes = path + "pre-analyses/gisaid_hcov-19.fasta"
    # new_genomes = path + "pre-analyses/new_genomes.fasta"
    # keep = path + 'config/keep.txt'
    # remove = path + "config/remove.txt"
    # outfile = path + "pre-analyses/temp_sequences.fasta"

    # inspect genome coverage
    genome_size = 29903
    max_gaps = 30
    min_size = genome_size - int(genome_size * (max_gaps / 100))

    # store only new sequences in a dictionary, ignoring existing ones
    print('\n### Loading new sequences, and reporting genome coverage\n')
    newly_sequenced = {}
    low_coverage = {}
    for fasta in SeqIO.parse(open(new_genomes), 'fasta'):
        id, seq = fasta.description, fasta.seq
        size = len(str(seq).replace('N', '').replace('-', ''))
        if size > min_size:
            coverage = str(round(size / genome_size, 3))
            print(id + ', coverage = ' + coverage + ' (PASS)')
            if id not in newly_sequenced.keys():  # avoid potential duplicates
                newly_sequenced[id] = str(seq)
        else:
            coverage = str(round(size / genome_size, 3))
            print(id + ', coverage = ' + coverage + ' (FAIL)')
            low_coverage[id] = coverage
    print('\nDone!\n')

    # create a list of sequences to be added in all instances
    keep_sequences = []
    for id in open(keep, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id not in newly_sequenced:
                if id not in keep_sequences:
                    keep_sequences.append(id)

    # create a list of sequences to be ignored in all instances
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)

    # export only sequences to be used in the nextstrain build
    c = 1
    print('\n### Exporting sequences\n')
    exported = []
    ignored = []
    all_sequences = []
    with open(outfile, 'w') as output:
        if genomes.split('.')[1] == 'json':
            print('\n### Exporting sequences\n')
            with open(genomes) as infile:
                c = 0
                for line in infile:
                    entry = json.loads(line)
                    id = entry['covv_virus_name']
                    id = id.split('|')[0].replace('hCoV-19/', '')
                    if len(id.split('/')) == 3:
                        country, index, year = id.split('/')
                    elif len(id.split('/')) == 4:
                        host, country, index, year = id.split('/')
                    country = country.replace(' ', '').replace('\'', '-').replace('_', '')
                    id = '/'.join([country, index, year])

                    all_sequences.append(id)
                    if id not in remove_sequences:
                        if id in keep_sequences:  # filter out unwanted sequences
                            seq = entry['sequence'].replace('\n', '')
                            if int(len(seq)) >= min_size:
                                entry = ">" + id + "\n" + seq.upper() + "\n"
                                exported.append(id)
                                output.write(entry)
                                print(str(c) + '. ' + id)
                                c += 1
                    else:
                        ignored.append(id)

        else:
            for fasta in SeqIO.parse(open(genomes), 'fasta'):
                id, seq = fasta.description, fasta.seq
                id = id.split('|')[0].replace('hCoV-19/', '')
                if len(id.split('/')) == 3:
                    country, index, year = id.split('/')
                elif len(id.split('/')) == 4:
                    host, country, index, year = id.split('/')
                country = country.replace(' ', '').replace('\'', '-').replace('_', '')
                id = '/'.join([country, index, year])
                all_sequences.append(id)

                if id not in remove_sequences:
                    if id in keep_sequences:  # filter out unwanted sequences
                        entry = ">" + id + "\n" + str(seq).upper() + "\n"
                        exported.append(id)
                        output.write(entry)
                        print(str(c) + '. ' + id)
                        c += 1
                else:
                    ignored.append(id)

        for id, seq in newly_sequenced.items():
            print('* ' + str(c) + '. ' + id)
            entry = ">" + id + "\n" + seq.upper() + "\n"
            exported.append(id)
            output.write(entry)
            c += 1
    print('\n- Done!\n')

    # mismatched sequence headers
    mismatch = [genome for genome in keep_sequences if genome not in all_sequences]
    if len(mismatch) + len(low_coverage) > 0:
        print('\n### WARNINGS!')

    if len(mismatch) > 0:
        print('\n### Possible sequence header mismatches\n')
        m = 1
        for id in mismatch:
            print(str(m) + '. ' + id)
            m += 1
    else:
        print('\nNo sequence name mismatches found...')

    if len(low_coverage) > 0:
        print('\n- Low quality sequences were ignored.\n')
        l = 1
        for id, coverage in low_coverage.items():
            print('\t' + str(l) + '. ' + id + ', coverage = ' + coverage + ' (FAIL)')
            l += 1

    print('\n\n\n### Final result\n')

    print('Lab file contains ' + str(len(newly_sequenced)) + ' high coverage sequences')
    print('Lab file contains ' + str(len(low_coverage)) + ' low coverage sequences, which were ignored')
    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(newly_sequenced)) + ' newly sequenced genomes were added')
    print(str(len(low_coverage)) + ' low coverage genomes were ignored')
    print(str(len(ignored)) + ' genomes were REMOVED according to remove.txt\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')
