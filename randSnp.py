import sys
import random
import string
prob_snp = 0.01
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


usage = ' randSnp [in.fastq] '
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print >>stderr, usage
        exit(1)
    #print "cmd line :   %s"%(usage)
    
    try:
        fp = open(sys.argv[1], 'r')
    except IOerror:
        print "can not open", sys.argv[1]
        exit(1)
    for name, seq, qual in readfq(fp):
        seqName = name.split('\t')[0]
        for i, nt in enumerate(seq):
            if random.random() <= prob_snp: #have a snp in pos i
                r = random.random()
                if r >0.0 and r < 0.333:#2 nucleoside type
                    snp_list = []
                    snp_list.append(string.upper(nt))
                    nucleoside = random.choice('ACGT')
                    while nucleoside == string.upper(nt):
                        nucleoside = random.choice('ACGT')
                    snp_list.append(nucleoside)
                    snp_list.sort()
                    snp_field = ''
                    #print snp_list
                    for nucleoside in snp_list:
                        snp_field += nucleoside +'/'
                    print "%s\t%u\t%s\t%s"%(seqName, i+1, seq[i], snp_field[:-1])
                elif r >0.333 and r < 0.666:#3 nucleoside type
                    snp_list = []
                    snp_list.append(string.upper(nt))
                    while len(snp_list) < 3:
                        nucleoside = random.choice('ACGT')
                        while nucleoside in snp_list: 
                            nucleoside = random.choice('ACGT')
                        snp_list.append(nucleoside)
                    snp_list.sort()
                    snp_field = ''
                    #print snp_list
                    for nucleoside in snp_list:
                        snp_field += nucleoside +'/'
                    print "%s\t%u\t%s\t%s"%(seqName, i+1, seq[i], snp_field[:-1])



                else:#4 nucleoside type
                    print "%s\t%u\t%s\t%s"%(seqName, i+1, seq[i], 'A/C/G/T')
    fp.close()    