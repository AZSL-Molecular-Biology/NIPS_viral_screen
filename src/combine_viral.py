#/bin/python

import argparse
import os
import sys
import json



def main(run: str, unit: str, sample: str, viral_file: str, factor: float, input_dir: str, output_file: str, verbose: bool=False):
    if verbose:
        print("Using:")
        print("\tViral File:\t{}".format(viral_file))
        print("\tFactor:\t{}".format(factor))
        print("\tInput directory:\t{}".format(input_dir))
        print("\tOutput file name:\t{}".format(output_file))

    # read viral file
    viral_dict = dict()
    family_dict = dict()
    with open(viral_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line.startswith("#"):
                line_list = line.split("\t")
                viral_dict[line_list[0]] = set()
                if len(line_list) > 2:
                    if line_list[2] not in family_dict:
                        family_dict[line_list[2]] = list()
                    family_dict[line_list[2]].append(line_list[0])
    if verbose:
        print(f"Found {len(viral_dict)} viral genomes")
        print(f"Found {len(family_dict)} families")

    # read all viral input files
    # store in read_dict: read names -> list of viral
    # store in viral_dict: viral -> list of read names 
    read_dict = dict()
    for viral_name in viral_dict:
        with open(f"{input_dir}/{viral_name}.highq_fragments.txt", "r") as f:
            for line in f:
                line = line.rstrip()
                line_split = line.split("\t")
                read_name = line_split[0]
                viral_dict[viral_name].add(read_name)
                if read_name not in read_dict:
                    read_dict[read_name] = list()
                read_dict[read_name].append(viral_name)

    if verbose:
        print("All viral files are parsed")

    # create counts for viral in viral_count_dict
    # viral_count_dict: viral -> count of uniq reads
    viral_count_dict = dict()
    for viral_name in viral_dict:
        viral_count_dict[viral_name] = 0
        for read in viral_dict[viral_name]:
            if len(read_dict[read]) == 1:
                viral_count_dict[viral_name] += 1

    # create counts for family
    family_read_dict = dict()
    for family in family_dict:
        family_read_dict[family] = 0
        read_set = set()
        for viral_name in family_dict[family]:
            for read in viral_dict[viral_name]:
                read_set.add(read)
        for read in read_set:
            for viral_name in read_dict[read]:
                matched_read = True
                if viral_name not in family_dict[family]:
                    matched_read = False
                if matched_read:
                    family_read_dict[family] += 1

    if verbose:
        print("All families are processed")

    # create json
    my_dict = dict()
    my_dict["run"] = run
    my_dict["unit"] = unit
    my_dict["sample"] = sample
    my_dict["normalization_factor"] = factor
    my_dict["viral_list"] = list()
    for viral_name in viral_dict:
        v_dict = dict()
        v_dict["name"] = viral_name
        v_dict["count"] = viral_count_dict[viral_name]
        v_dict["normalized"] = viral_count_dict[viral_name]/factor
        my_dict["viral_list"].append(v_dict)
    my_dict["family_list"] = list()
    for family in family_read_dict:
        f_dict = dict()
        f_dict["name"] = family
        f_dict["contains"] = family_dict[family]
        f_dict["count"] = family_read_dict[family]
        f_dict["normalized"] = family_read_dict[family]/factor
        my_dict["family_list"].append(f_dict)

    with open(output_file, 'w') as f:
        f.write(json.dumps(my_dict, indent=4))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The AZ Delta viral combiner')
    parser.add_argument("run", help='The run')
    parser.add_argument("unit", help='The unit')
    parser.add_argument("sample", help='The sample')
    parser.add_argument("viral_file", help='The file with viral information')
    parser.add_argument("factor", help='The factor of normalization', type=float)
    parser.add_argument("input_dir", help='The input of all files of viral sequences')
    parser.add_argument("output_file", help='The name of output json')
    parser.add_argument('-v', "--verbose", help='verbose', action='store_true')
    
    args = parser.parse_args()
    run = args.run
    unit = args.unit
    sample = args.sample
    viral_file = args.viral_file
    factor = args.factor
    input_dir = args.input_dir
    output_file = args.output_file
    verbose = args.verbose

    main(run, unit, sample, viral_file, factor, input_dir, output_file, verbose=verbose)