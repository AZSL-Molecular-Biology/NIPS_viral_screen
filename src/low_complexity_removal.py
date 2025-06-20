import gzip
import argparse
import math

class Read():

    def __init__(self, name, read1, qual1, read2, qual2):
        self._name = name
        self._read1 = read1
        self._read2 = read2
        self._qual1 = qual1
        self._qual2 = qual2

    def get_qual(self) -> float:
        qual = 0
        for i in self._qual1:
            qual += (ord(i) - 33)
        for i in self._qual2:
            qual += (ord(i) - 33)
        return qual / (len(self._qual1 + self._qual2))

    def _get_dust(self, sequence, kmer: int=3) -> int:
        nuc_freq = dict()

        for i in range(0, len(sequence) - kmer):
            sub = sequence[i:i+kmer]
            nuc_freq[sub] = sequence.count(sub)

        dust = 0
        for nuc in nuc_freq:
            freq = nuc_freq[nuc]
            dust = dust + (freq * (freq-1) * 0.5)
        
        return dust * 1/(len(sequence) - (kmer + 1))

    def get_max_dust(self):
        d1 = self._get_dust(self._read1)
        d2 = self._get_dust(self._read2)
        return max([d1, d2])

    def __lt__(self, o):
        return self.get_qual() < o.get_qual()

def main(fastq1: str, fastq2: str, outfastq1, outfastq2, verbose: bool=False):
    read_list = list()
    with gzip.open(fastq1,'rt') as f1:
        with gzip.open(fastq2, 'rt') as f2:
            while True:
                line = f1.readline()
                if not line:
                    break
                name = line.strip()
                f2.readline()
                read1 = f1.readline().strip()
                read2 = f2.readline().strip()
                f1.readline()
                f2.readline()
                qual1 = f1.readline().strip()
                qual2 = f2.readline().strip()
                read_list.append(Read(name, read1, qual1, read2, qual2))
    
    if verbose:
        print("Found {} reads".format(len(read_list)))

    # find duplicates
    dup_dict = dict()
    for read in read_list:
        dup_read = "{}-{}".format(read._read1, read._read2)
        if dup_read not in dup_dict:
            dup_dict[dup_read] = list()
        dup_dict[dup_read].append(read)
    if verbose:
        print("Found {} uniques".format(len(dup_dict)))
        has_multiple_counts = 0
        most = 0
        for r in dup_dict:
            if len(dup_dict[r]) > 1:
                has_multiple_counts += 1
                if len(dup_dict[r]) > most:
                    most = len(dup_dict[r])
        print("Found {} reads with multiple occurances, most was {}".format(has_multiple_counts, most))

    # keep best reads (best duplicates)
    dedup_list = list()
    for r in dup_dict:
        dedup_list.append(sorted(dup_dict[r])[-1])

    if verbose: 
        print("Kept {} reads to remove low complexity".format(len(dedup_list)))

    dust_list = list()
    for read in dedup_list:
        if read.get_max_dust() < 1:
            dust_list.append(read)
            
    if verbose:
        print("Kept {} reads after DUST filtering".format(len(dust_list)))

    with gzip.open(outfastq1,'wt') as f1:
        with gzip.open(outfastq2, 'wt') as f2:
            for read in dust_list:
                f1.write(read._name + '\n')
                f2.write(read._name + '\n')
                f1.write(read._read1 + '\n')
                f2.write(read._read2 + '\n')
                f1.write("+" + '\n')
                f2.write("+" + '\n')
                f1.write(read._qual1 + '\n')
                f2.write(read._qual2 + '\n')

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The AZ Delta viral fastq filter')
    parser.add_argument("fastq1", help='The first fastq file')
    parser.add_argument("fastq2", help='The second fastq file')
    parser.add_argument("outfastq1", help='The first output fastq file')
    parser.add_argument("outfastq2", help='The second output fastq file')
    parser.add_argument('-v', "--verbose", help='verbose', action='store_true')

    args = parser.parse_args()
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    outfastq1 = args.outfastq1
    outfastq2 = args.outfastq2
    verbose = args.verbose

    main(fastq1, fastq2, outfastq1, outfastq2, verbose=verbose)