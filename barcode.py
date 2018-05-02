#!/sas5/home/jiangdawei/python3/bin/python3
import random
import time
from functools import wraps
import argparse
import copy
import gzip
import numpy as np
from collections import namedtuple
Match = namedtuple('Match', 'a, b, size')

'''用于测试函数运行时间'''

def fn_timer(function):
	@wraps(function)
	def function_timer(*args, **kwargs):
		t0 = time.time()
		result = function(*args, **kwargs)
		t1 = time.time()
		print("Total time running %s: %s seconds" % (function.__name__, str(t1-t0)))
		return result
	return function_timer

'''
The class SequenceMatcher is modified from difflib model of python.
An algorithm published in the late 1980’s by Ratcliff and Obershelp 
under the hyperbolic name “gestalt pattern matching.” 
'''
#@fn_timer
class SequenceMatcher:
	def __init__(self, a='', b=''):
		self.a =self.b = None
		self.set_seqs(a,b)

	def set_seqs(self, a, b):
		self.set_seq1(a)
		self.set_seq2(b)

	def set_seq1(self, a):
		if a is self.a:
			return
		self.a = a
		self.matching_blocks = self.opcodes = None

	def set_seq2(self, b):
		if b is self.b:
			return
		self.b = b
		self.matching_blocks = self.opcodes = None
		self.fillbcount = None
		self.__chain_b()

	def __chain_b(self):
		b = self.b
		self.b2j = b2j = {}

		for i, elt in enumerate(b):
			indices = b2j.setdefault(elt, [])
			indices.append(i)

	def find_longest_match(self, alo, ahi, blo, bhi):
		a, b, b2j = self.a ,self.b, self.b2j
		besti, bestj, bestsize, = alo, blo, 0
		j2len = {}
		nothing = []
		for i in range(alo, ahi):
			j2lenget = j2len.get
			newj2len = {}
			for j in b2j.get(a[i], nothing):
				if j < blo:
					continue
				if j >=bhi:
					break
				k = newj2len[j] = j2lenget(j-1,0) +1
				if k > bestsize:
					besti, bestj, bestsize = i-k+1, j-k+1, k
			j2len = newj2len
		return Match(besti, bestj, bestsize)
	
	def get_matching_blocks(self):
		if self.matching_blocks is not None:
			return self.matching_blocks
		la, lb = len(self.a), len(self.b)
		queue = [(0, la, 0, lb)]
		matching_blocks = []
		while queue:
			alo, ahi, blo, bhi = queue.pop()
			i, j, k = x = self.find_longest_match(alo, ahi, blo, bhi)
			if k:
				matching_blocks.append(x)
				if alo < i and blo < j:
					queue.append((alo, i, blo, j))
				if i+k < ahi and j+k < bhi:
					queue.append((i+k, ahi, j+k, bhi))
		matching_blocks.sort()
		
		i1 = j1 = k1 =0
		non_adjacent = []
		for i2, j2, k2 in matching_blocks:
			if i1 + k1 == i2 and j1+ k1 ==j2:
				k1 += k2
			else:
				if k1:
					non_adjacent.append((i1, j1, k1))
				i1, j1, k1 = i2, j2 ,k2
		if k1:
			non_adjacent.append((i1, j1, k1))

		non_adjacent.append((la, lb, 0))
		self.matching_blocks = list(map(Match._make, non_adjacent))
		return self.matching_blocks

	def get_opcodes(self):
		if self.opcodes is not None:
			return self.opcodes
		i = j = 0
		self.opcodes = answer = []
		for  ai, bj, size in self.get_matching_blocks():
			tag = ''
			if i < ai and j < bj:
				tag = 'replace'
			elif i < ai:
				tag = 'delete'
			elif j < bj:
				tag = 'insert'
			if tag:
				answer.append((tag, i, ai, j, bj))
			i, j = ai+size, bj+size
			if size:
				answer.append(('equal', ai, i, bj, j))
		return answer

	#@fn_timer
	def get_matchs(self):
		matches = sum(triple[-1] for triple in self.get_matching_blocks())
		return matches
	
	def ratio(self):
		matches = sum(triple[-1] for triple in self.get_matching_blocks())
		length = len(self.b) * 2
		if length:
			return 2.0 * matches / length
		return 1.0

class sample_info():
	def __init__(self, name, ratio, r1_ratio, r2_ratio, barcode1, barcode2):
		self.name = name
		self.ratio = ratio
		self.r1_ratio = r1_ratio
		self.r2_ratio = r2_ratio
		self.barcode1 = barcode1
		self.barcode2 = barcode2
	def __repr__(self):
		return repr((self.name, self.ratio, self.r1_ratio, self.r2_ratio, self.barcode1, self.barcode2))
#@fn_timer
class PrimerFinder:
	
	def __init__(self, r1, r2, pf, pr):
		"""Construct a PrimerFinder.
		Arg r1, r2 are sequence of pair reads.r1 and r2 are not ordered.
		Arg pf, pr are primer sequence.
		"""
		self.r1 = r1
		self.r2 = r2
		self.pf = pf
		self.pr = pr
	
	def find_position(self, SequenceMatcher_result, primer_len):
		for maptype, r_start, r_end, p_start, p_end in SequenceMatcher_result.get_opcodes():
			if p_start == 0 and p_end > p_start:
				match_start = r_start
			if p_end == primer_len and p_end > p_start:
				match_end = r_end
		return np.array([match_start, match_end])


	def primer_match(self):
		forward_1 = SequenceMatcher(self.r1[6:43], self.pf)
		forward_2 = SequenceMatcher(self.r2[6:55], self.pr)
		reverse_1 = SequenceMatcher(self.r1[6:55], self.pr)
		#print(self.r1[6:55], self.pr)
		#print(reverse_1.get_opcodes())
		reverse_2 = SequenceMatcher(self.r2[6:43], self.pf)
		#print(self.r2[6:43], self.pf)
		#print(reverse_2.get_opcodes())
		forward_matchs = forward_1.get_matchs() + forward_2.get_matchs()
		reverse_matchs = reverse_1.get_matchs() + reverse_2.get_matchs()
		#print(reverse_1.get_matchs())
		#print(reverse_2.get_matchs())
		if(forward_matchs > reverse_matchs):
			orientation = "forward"
			r1_pos = self.find_position(forward_1, len_primerf)
			r1_pos = r1_pos + 6
			r2_pos = self.find_position(forward_2, len_primerr)
			r2_pos = r2_pos + 6
		else:
			orientation = "reverse"
			r1_pos = self.find_position(reverse_1, len_primerr)
			#print(r1_pos)
			r1_pos = r1_pos + 6
			r2_pos = self.find_position(reverse_2, len_primerf)
			#print(r2_pos)
			r2_pos = r2_pos + 6
		#print(orientation, r1_pos, r2_pos)
		return(orientation, r1_pos, r2_pos)


#@fn_timer
class BarcodeFinder:
	
	def __init__(self, orientation, seq1_bar, seq2_bar, samplebar):
		"""Construct a BarcodeFinder.
		Arg seq1_bar, seq2_bar are sequence of pair reads before primers.
		Arg sample is tuple storing barcode+spacer+primer sequence of samples.
		"""
		self.orientation = orientation
		self.seq1_bar = seq1_bar
		self.seq2_bar = seq2_bar
		self.samplebar = samplebar
	
	#@fn_timer
	def sample_ratio(self,sample):
		"""Calculate match ratio of a sample"""
		r1 = SequenceMatcher(self.seq1_bar, sample.barcode1)
		r2 = SequenceMatcher(self.seq2_bar, sample.barcode2)
		#print(self.seq1_bar, sample.barcode1, r1.get_matchs(), self.seq2_bar, sample.barcode2,r2.get_matchs())
		sample.r1_ratio = r1.get_matchs()/(len(self.seq1_bar)+len(sample.barcode1)) * 2
		sample.r2_ratio = r2.get_matchs()/(len(self.seq2_bar)+len(sample.barcode2)) * 2
		sample.ratio = (sample.r1_ratio + sample.r2_ratio) / 2
		return sample
	#@fn_timer
	def samplematch(self):
		"""Find maximum match ratio of all sample"""
		samplelist = []
		for k, s in self.samplebar.items():
			f, r = k.split("-")
			if self.orientation == "forward":
				sample = sample_info(s, 0, 0, 0, f, r)
			else:
				sample = sample_info(s, 0, 0, 0, r, f)
			sample = self.sample_ratio(sample)
			if sample.r1_ratio == 1 and sample.r2_ratio == 1:
				#print (sample, self.seq1_bar, self.seq2_bar)
				return [sample]
			else:
				samplelist.append(sample)
		newlist = sorted(samplelist, key=lambda sample:sample.ratio)
		#print(newlist[-1], self.seq1_bar, self.seq2_bar)
		return newlist

#@fn_timer
class qc:
	def __init__(self, seq1, seq2, q1, q2, qc_info, low_quality_threshold):
		self.seq1 = seq1
		self.seq2 = seq2
		self.q1 = q1
		self.q2 = q2
		self.qc_info = qc_info
		self.low_quality_threshold = low_quality_threshold
	
	#@fn_timer
	def ascii_to_phred33(self, q):
		return np.fromstring(q,dtype=np.int8) - 33
	
	#@fn_timer
	def qual_stat(self, qual_threshold, qual):
		qual = self.ascii_to_phred33(qual)
		compare = qual >= qual_threshold
		qual_num = np.sum(compare==True)
		return qual_num
	
	#@fn_timer
	def trim(self, qual_threshold, qual):
		qual = self.ascii_to_phred33(qual)
		compare = qual >= qual_threshold
		idx, = compare.nonzero()
		if len(idx > 0 ):
			trim_pos = idx[-1] + 1
			return trim_pos
		else:
			return 0
	
	#@fn_timer
	def gc_stat(self, sequence):
		g_num = sequence.count("G")
		c_num = sequence.count("C")
		total_num = len(sequence)
		return g_num+c_num, total_num

	#@fn_timer
	def qc_stat(self):
		self.qc_info["read1"]["q20_base"] += self.qual_stat(20, self.q1)
		self.qc_info["read2"]["q20_base"] += self.qual_stat(20, self.q2)
		self.qc_info["read1"]["q30_base"] += self.qual_stat(30, self.q1)
		self.qc_info["read2"]["q30_base"] += self.qual_stat(30, self.q2)
		gc, total = self.gc_stat(self.seq1)
		self.qc_info["read1"]["total_base"] += total
		self.qc_info["read1"]["gc_base"] += gc
		gc, total = self.gc_stat(self.seq2)
		self.qc_info["read2"]["total_base"] += total
		self.qc_info["read2"]["gc_base"] += gc
		read1_trim_pos = self.trim(self.low_quality_threshold, self.q1)
		read2_trim_pos = self.trim(self.low_quality_threshold, self.q2)
		q = self.q1 + self.q2
		high_percent = (self.qual_stat(self.low_quality_threshold, q))/len(q)
		return self.qc_info, read1_trim_pos, read2_trim_pos, high_percent

#@fn_timer
def readfq(f1, f2):
	while True:
		name1 = f1.readline().decode('UTF-8').strip()
		name2 = f2.readline().decode('UTF-8').strip()
		seq1 = f1.readline().decode('UTF-8').strip()
		seq2 = f2.readline().decode('UTF-8').strip()
		info1 = f1.readline().decode('UTF-8').strip()
		info2 = f2.readline().decode('UTF-8').strip()
		q1 = f1.readline().decode('UTF-8').strip()
		q2 = f2.readline().decode('UTF-8').strip()
		if not q1 or not q2:
			break
		yield (name1, name2, seq1, seq2, info1, info2, q1, q2)

#@fn_timer
def readsample(barcodefile, primerfile, samplefile):
	with open(barcodefile) as b, open(primerfile) as p, open(samplefile) as sam:
		barcodedict = {}
		#sample = []
		sample = {}
		primerf, primerr = p.readline().strip().split('\t')
		global len_primerf
		global len_primerr
		len_primerf = len(primerf)
		len_primerr = len(primerr)
		for line in b.readlines():
			t, barcode, spacer = line.strip().split('\t')
			if t.startswith('#'):
				continue
			else:
				barcodedict[barcode] = spacer
				'''
				此处不考虑匹配pcr引物
				if t.startswith('SF'):
					barcodedict[barcode] = spacer + primerf
				elif t.startswith('SR'):
					barcodedict[barcode] = spacer + primerf
				'''
		for line in sam.readlines():
			s, b1, b2 = line.strip().split('\t')
			#sample.append((s, barcodedict[b1], barcodedict[b2]))
			key = barcodedict[b1]+"-"+barcodedict[b2]
			sample[key] = s
		#print(sample)
	return  primerf, primerr, sample

#def _contiguous_regions()

def GZOut(outfile1, outfile2, name1, name2, seq1, seq2, info1, info2, q1, q2):
	with gzip.open(outfile1,'a') as o1, gzip.open(outfile2,'a') as o2:
					line1 = name1+'\n'+seq1+'\n'+info1+'\n'+q1+'\n'
					line2 = name2+'\n'+seq2+'\n'+info2+'\n'+q2+'\n'
					o1.write(line1.encode('UTF-8'))
					o2.write(line2.encode('UTF-8'))

def Out(outfile1, outfile2, name1, name2, seq1, seq2, info1, info2, q1, q2):
	with open(outfile1,'a') as o1, open(outfile2,'a') as o2:
					line1 = name1+'\n'+seq1+'\n'+info1+'\n'+q1+'\n'
					line2 = name2+'\n'+seq2+'\n'+info2+'\n'+q2+'\n'
					o1.write(line1)
					o2.write(line2)
def main():
	
	parser = argparse.ArgumentParser(description='Split sample by barcodes from 16s sequence data.')
	parser.add_argument('--fq1', help='The first input fastq file. Compressed format - gz')
	parser.add_argument('--fq2', help='The second input fastq file. Compressed format - gz')
	parser.add_argument('--barcodefile', help='Bardode and Spacer sequence')
	parser.add_argument('--primerfile', help='Primer sequence')
	parser.add_argument('--samplefile', help='Sample information')
	parser.add_argument('--minratio', help='Min match ratio (Average value for two barcode).')
	parser.add_argument('--low_quality_threshold', help='Low quality threshold.')
	parser.add_argument('--high_quality_percent', help='Minimum High quality base percent of reads pair.')
	parser.add_argument('--outpath', help='Output path ')
	args = parser.parse_args()
	fq1 = args.fq1
	fq2 = args.fq2
	barcodefile = args.barcodefile
	primerfile = args.primerfile
	samplefile = args.samplefile
	minratio = float(args.minratio)
	low_quality_threshold = int(args.low_quality_threshold)
	high_quality_percent = float(args.high_quality_percent)
	outpath =args.outpath + '/'

	'''Read barcodefile, primerfile, samplefile.'''
	primerf, primerr, samplebar = readsample(barcodefile, primerfile, samplefile)
	#print (primerf, primerr, samplebar)
	
	''' Set Initial Value. '''
	i = 0
	total = no_match = multiple_match = quality_filter_num = 0
	readsnum = {}
	qc_info = {
			"read1":{
				"q20_base":0,
				"q30_base":0,
				"gc_base":0,
				"total_base":0,
				},
			"read2":{
				"q20_base":0,
				"q30_base":0,
				"gc_base":0,
				"total_base":0,
				}
	}

	'''Open fq.gz file'''
	f1 = gzip.open (fq1)
	f2 = gzip.open (fq2)
	read1_trim_pos, read2_trim_pos = 0, 0
	for name1, name2 ,rawseq1, rawseq2, info1, info2, rawq1, rawq2 in readfq(f1, f2):
		i +=1
		qc_info, read1_trim_pos, read2_trim_pos, high_percent = qc(rawseq1, rawseq2, rawq1, rawq2, qc_info, low_quality_threshold).qc_stat()
		if high_percent < high_quality_percent:
			quality_filter_num +=1
			#outfile1 = outpath + 'low_quality_R1.fq.gz'
			#outfile2 = outpath + 'low_quality_R2.fq.gz'
			outfile1 = outpath + 'low_quality_R1.fq'
			outfile2 = outpath + 'low_quality_R2.fq'
			seq1 = rawseq1
			seq2 = rawseq2
			q1 = rawq1
			q2 = rawq2
		elif name1.startswith('@'):
			orientation, r1_pos, r2_pos = PrimerFinder(rawseq1, rawseq2, primerf, primerr).primer_match()
			seq1_bar = rawseq1[0:r1_pos[0]]
			seq2_bar = rawseq2[0:r2_pos[0]]
			seq1 = rawseq1[r1_pos[1]:read1_trim_pos]
			seq2 = rawseq2[r2_pos[1]:read2_trim_pos]
			q1 = rawq1[r1_pos[1]:read1_trim_pos]
			q2 = rawq2[r2_pos[1]:read2_trim_pos]
			k = seq1_bar+"-"+seq2_bar
			#print(k)
			if k in samplebar:
				best_name = samplebar[k]
				readsnum[best_name] = readsnum.get(best_name, 0)
				readsnum[best_name] += 1
				#outfile1 = outpath + best_name+'_R1.fq.gz'
				#outfile2 = outpath + best_name+'_R2.fq.gz'
				outfile1 = outpath + best_name+'_R1.fq'
				outfile2 = outpath + best_name+'_R2.fq'
			else:
				sample_match = BarcodeFinder(orientation, seq1_bar, seq2_bar, samplebar).samplematch()
				if sample_match[-1].ratio == 1:
					best_match = sample_match[-1]
					best_name=best_match.name
					readsnum[best_name] = readsnum.get(best_name, 0)
					readsnum[best_name] += 1
					#outfile1 = outpath + best_name+'_R1.fq.gz'
					#outfile2 = outpath + best_name+'_R2.fq.gz'
					outfile1 = outpath + best_name+'_R1.fq'
					outfile2 = outpath + best_name+'_R2.fq'
				elif sample_match[-1].r1_ratio > minratio and sample_match[-1].r2_ratio > minratio:
					if sample_match[-1].ratio > sample_match[-2].ratio:
						best_match = sample_match[-1]
						best_name=best_match.name
						readsnum[best_name] = readsnum.get(best_name, 0)
						readsnum[best_name] += 1
						#outfile1 = outpath + best_name+'_R1.fq.gz'
						#outfile2 = outpath + best_name+'_R2.fq.gz'
						outfile1 = outpath + best_name+'_R1.fq'
						outfile2 = outpath + best_name+'_R2.fq'
					elif sample_match[-1].ratio == sample_match[-2].ratio:
						multiple_match += 1
						#outfile1 = outpath + 'multiple_match_R1.fq.gz'
						#outfile2 = outpath + 'multiple_match_R2.fq.gz'
						outfile1 = outpath + 'multiple_match_R1.fq'
						outfile2 = outpath + 'multiple_match_R2.fq'
						seq1 = rawseq1
						seq2 = rawseq2
						q1 = rawq1
						q2 = rawq2
				else:
					no_match += 1
					#outfile1 = outpath + 'no_match_R1.fq.gz'
					#outfile2 = outpath + 'no_match_R2.fq.gz'
					outfile1 = outpath + 'no_match_R1.fq'
					outfile2 = outpath + 'no_match_R2.fq'
					seq1 = rawseq1
					seq2 = rawseq2
					q1 = rawq1
					q2 = rawq2
		#GZOut(outfile1, outfile2, name1, name2, seq1, seq2, info1, info2, q1, q2)
		Out(outfile1, outfile2, name1, name2, seq1, seq2, info1, info2, q1, q2)
	f1.close()
	f2.close()
	total = i
	
	qcfile = outpath + 'qc.xls'
	with open(qcfile,'w') as o:
		clean_num = total - quality_filter_num
		o.write('Total_read_pair_num\tQuality_filter_num\tPercent\tClean_num\tPercent\tRead1_Base\t%GC\t%Q20\t%Q30\tRead2_Base\t%GC\t%Q20\t%Q30\tTotal_Base\t%GC\t%Q20\t%Q30\n')
		o.write('%d\t' % (total))
		o.write('%d\t%.3f\t' % (quality_filter_num, quality_filter_num/total*100))
		o.write('%d\t%.3f\t' % (clean_num, clean_num/total*100))
		read1 = qc_info['read1']
		TotalBase1 = read1['total_base']
		o.write('%d\t%.3f\t%.3f\t%.3f\t' % (TotalBase1, read1['gc_base']/TotalBase1*100, read1['q20_base']/TotalBase1*100, read1['q30_base']/TotalBase1*100))
		read2 = qc_info['read2']
		TotalBase2 = read2['total_base']
		o.write('%d\t%.3f\t%.3f\t%.3f\t' % (TotalBase2, read2['gc_base']/TotalBase2*100, read2['q20_base']/TotalBase2*100, read2['q30_base']/TotalBase2*100))
		TotalBase = TotalBase1 + TotalBase2
		o.write('%d\t%.3f\t%.3f\t%.3f\n' % (TotalBase, (read1['gc_base']+read2['gc_base'])/TotalBase*100, (read1['q20_base']+read2['q20_base'])/TotalBase*100, (read1['q30_base']+read2['q30_base'])/TotalBase*100))

	outfile = outpath + 'split.xls'
	with open(outfile,'w') as o:
		split_num = clean_num - multiple_match - no_match
		o.write('Clean_num\t%d\n' % (clean_num))
		o.write('Splitted_read_pair_num\t%d\t%.3f\n' % (split_num, split_num/clean_num*100))
		o.write('Multiple_match\t%d\t%.3f\n' % (multiple_match, multiple_match/clean_num*100))
		o.write('No_match\t%d\t%.3f\n' % (no_match, no_match/clean_num*100))
		for j,j_num in sorted(readsnum.items(),key = lambda x:x[1],reverse = True):
			o.write(j+'\t%d\t%.3f\n' % (j_num, j_num/clean_num*100))

if __name__ == '__main__':
	main()
