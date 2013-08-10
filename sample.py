def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

import sys
import random
N = 2*1024*1024 #2M
usage = 'sample.py [in1.fastq] [in2.fastq]'
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print usage
        exit(1)
    #count reads number in in1.fastq
    fp = open(sys.argv[1], 'r')
    n = 0
    for name, seq, qual in readfq(fp):
        n += 1
    fp.close()
    if n < N:
        print >>stderr, 'seq number smaller than sample number!'
        exit(1)
    #sample
    population = [i for i in range(n)]
    sample = random.sample(population, N)
    #choice reads from in1.fastq
    fp = open(sys.argv[1], 'r')
    fp_o = open(sys.argv[1]+'.sample', 'w')
    for i, (name, seq, qual) in enumerate(readfq(fp)):
        if i in sample:
            print >>fp_o, "@"+name
            print >>fp_o, seq
            if qual != None:
                print >>fp_o, '+'
                print >>fp_o, qual
    fp_o.close()
    fp.close()
    #choice reads from in2.fastq
    fp = open(sys.argv[2], 'r')
    fp_o = open(sys.argv[2]+'.sample', 'w')

    for i, (name, seq, qual) in enumerate(readfq(fp)):
        if i in sample:
            print >>fp_o, "@"+name
            print >>fp_o, seq
            if qual != None:
                print >>fp_o, '+'
                print >>fp_o, qual
    fp_o.close()
    fp.close()
