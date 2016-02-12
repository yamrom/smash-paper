import pysam
import numpy as np
import sys


## proceeds through a samfile
## aggregating alignments with the same qname
## and organizing by read1 and read2
## supports python iteration methods
class readIterator():
    def __init__(self, infilename):
        self.sam = pysam.Samfile(infile, 'rb')
        self.current     = self.sam.next()
        self.currentName = self.current.qname
        
    ## get the next entry
    def next(self):
        if self.currentName == None:
            raise StopIteration
        ans1, ans2 = [], []
        while self.current.qname == self.currentName:
            if self.current.is_read1:
                ans1.append(self.current)
            else:
                ans2.append(self.current)
            
            try:
                self.current = self.sam.next()
            except:
                self.current = None
                self.currentName = None

                return ans1, ans2
        self.currentName = self.current.qname
        return ans1, ans2
    
    ## return a reference name from a tid 
    ## (integer mapping of chrom names)
    def getrname(self, tid):
        if tid == -1:
            return "*"
        else:
            return self.sam.getrname(tid)
                
    def __iter__(self):
        self.sam.reset()
        self.current     = self.sam.next()
        self.currentName = self.current.qname
        return self
    
    ## close the samfile
    def close(self):
        self.sam.close()

## given a matchCode from a read set 
## (integer count of the hits overlapping those bases in the read)
## return the proportion of unique bases out of total matched bases in the read.
def getRatio(x, matchCode):
    mcount = x.qlen
    if x.is_reverse:
        end   = x.rlen - x.qstart
        start = x.rlen - x.qend
    else:
        start = x.qstart
        end   = x.qend
    oneCount = np.sum(matchCode[start:end] == 1)
    return oneCount / float(mcount)

## return all unique base ratios for a read set
def getRatios(readSet):
    matchCode = matchCounter(readSet)
    return [getRatio(x, matchCode) for x in readSet]
    

## remove all hits that are unmapped or below the threshold for minMatch     
def matchCountFilter(readSet, minMatch):
    ans = []
    for read in readSet: 
        if (not read.is_unmapped)  and (read.qlen >= minMatch):
            ans.append(read)
    return ans

## remove hits without sufficient excess mappability
def excessMappabilityFilter(readSet, minExcessMappability):
    ans = []
    for read in readSet: 
        if not read.is_unmapped:
            maxMappability = max(read.opt("L0"), read.opt("R0"))
            # print "%d maxMappability\t%d qlen" % (maxMappability, read.qlen)
            if (read.qlen - maxMappability >= minExcessMappability):
                ans.append(read)
    return ans


## for each read in the readSet
def matchCounter(readSet):
    if len(readSet) == 0:
        return []
    else:
        rlen = readSet[0].rlen
        ans = np.zeros(rlen, dtype=int)
        for x in readSet:
            if x.is_reverse:
                end   = rlen - x.qstart
                start = rlen - x.qend
            else:
                start = x.qstart
                end   = x.qend
            ans[start:end] += 1
        return ans

    
def matchCountFilterPair(readSet1, readSet2, minMatch):
    return matchCountFilter(readSet1, minMatch), matchCountFilter(readSet2, minMatch)


def excessMappabilityFilterPair(readSet1, readSet2, minExcessMappability):
    return excessMappabilityFilter(readSet1, minExcessMappability), excessMappabilityFilter(readSet2, minExcessMappability)



def getKey(r1_chrom, r1_pos, r1_order, r2_chrom, r2_pos, r2_order):
    chrom = []
    pos = []
    for x in r1_order:
        chrom.append(r1_chrom[x])
        pos.append(r1_pos[x])
    for x in r2_order:
        chrom.append(r2_chrom[x])
        pos.append(r2_pos[x])
    return tuple(chrom), tuple(pos)


#infile = "/data/safe/levy/smash/12871.s/longmem.ns.bam"

headings = ["read_id", 'read_index', 'hit_index', 'chrom', 'pos', 'reverse', 'read_len', 'hit_offset', 'match_len', 'umatch', 'excess'] 

infile = sys.argv[1]
minMatch = int(sys.argv[2])
minRatio = float(sys.argv[3])
hitWindow = int(sys.argv[4])
minExcessMappability = int(sys.argv[5])
iterator = readIterator(infile)

## print out the headings
print "\t".join(headings)

## keep a set of keys to watch for duplicates
dupeSet = set()
n_dupe = 0;
n_non_dupe = 0;

## for each read name in the datafile
for reads1, reads2 in iterator:
    ## get name
    readID = reads1[0].qname
    reads1, reads2 = excessMappabilityFilterPair(reads1, reads2, minExcessMappability)

    ## and filter lists for minimum match length    
    reads1, reads2 = matchCountFilterPair(reads1, reads2, minMatch)
    ## if the resulting lists are not empty
    if (len(reads1) > 0) or (len(reads2) > 0):
        ## get ratios for 1 and 2
        ratio1 = getRatios(reads1)
        ratio2 = getRatios(reads2)
                
        read1_info = []     ## store the mapping info
        read2_info = []     
        r1_chrom = []       ## chrom
        r1_pos   = []       
        r2_chrom = []       ## position
        r2_pos   = []
        r1_hitIndex = []    ## and hit index
        r2_hitIndex = []
        
        ## for each read in read1
        for x, ratio in zip(reads1, ratio1):
            ## if the ratio is good
            if ratio >= minRatio:
                ## get and save info
                hitIndex = x.opt("HI")
                info = (int(x.is_read2)+1, hitIndex, iterator.getrname(x.tid), x.pos, int(x.is_reverse), x.rlen, x.qstart, x.qlen, int(np.round(x.qlen*ratio)), x.qlen - max(x.opt("L0"), x.opt("R0")))
                r1_chrom.append(x.tid)
                r1_pos.append(x.pos)
                r1_hitIndex.append(hitIndex)                        
                read1_info.append(info)

        ## convert chrom and pos into numpy array for testing overlaps in read2
        r1_chrom = np.array(r1_chrom)
        r1_pos   = np.array(r1_pos) 
        
        ## for each read in read2
        for x, ratio in zip(reads2, ratio2):
            if ratio >= minRatio:
                position = x.pos
                # check if it is close to anything from read1
                matchRead1 = np.sum(np.logical_and(r1_chrom == x.tid,
                                                   np.abs(r1_pos - x.pos) < hitWindow))
                ## if not, save info                
                if matchRead1 == 0:                    
                    hitIndex = x.opt("HI")
                    info = (int(x.is_read2)+1, hitIndex, iterator.getrname(x.tid), x.pos, int(x.is_reverse), x.rlen, x.qstart, x.qlen, int(np.round(x.qlen*ratio)), x.qlen - max(x.opt("L0"), x.opt("R0")))

                    r2_hitIndex.append(hitIndex)
                    r2_chrom.append(x.tid)
                    r2_pos.append(x.pos)
                    
                    read2_info.append(info)
        
        r2_chrom = np.array(r2_chrom)
        r2_pos   = np.array(r2_pos)
        
        ## order by hit index
        r1_order = np.argsort(r1_hitIndex)
        r2_order = np.argsort(r2_hitIndex)
        ## make a key for dupeSet
        key = getKey(r1_chrom, r1_pos, r1_order, r2_chrom, r2_pos, r2_order)
        if not key in dupeSet:
            ## if not in the set, add the key and print output 
            dupeSet.add(key)
            for index in r1_order:
                print "%s\t%s" % (readID, "\t".join(map(str, read1_info[index])))
            
            for index in r2_order:
                print "%s\t%s" % (readID, "\t".join(map(str, read2_info[index])))
            n_non_dupe += 1
        else:
            n_dupe += 1

print "%d dupes\t%d non-dupes" % (n_dupe, n_non_dupe)

iterator.close()
