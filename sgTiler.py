#!/usr/bin/python

import sys, argparse
import math
import os
import subprocess
import shlex
import threading
import logging as log
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='The tool will identify tiling sgRNAs in a sequence', add_help=True)
parser.add_argument('-i', '--input', help = "Input fasta file", dest="FASTA")
parser.add_argument('--output', default='out', help="Output prefix")
parser.add_argument('--pam', dest='PAM', default='NGG', help='PAM sequence. Default: NGG')
parser.add_argument('--gtf', dest='GTF', help='Path to annotation gtf file')
parser.add_argument('--gc-min', type=int, default=20, help="Minimum percet GC content")
parser.add_argument('--gc-max', type=int, default=80, help="Maximum percent GC content")
parser.add_argument('-n', '--nthreads', type=int, default=4, help='No. of threads. Default: 4')
parser.add_argument('--handle', default='GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCG', help='dCas9 handle sequence. Default: GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCG')
parser.add_argument('--terminator', default='UUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU', help='Terminator sequence. Default: UUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU')
parser.add_argument('--strand', default='both', help='Strand to find the sgRNAs. Options: positive, negative or both. Default: both')
parser.add_argument('--length', '-l', type=int, default=19, help='Length of sgRNA without the PAM. Default: 19')
parser.add_argument('--bowtie-index', help='Bowtie index')
parser.add_argument('--dhs', help='Path to DHS file')
parser.add_argument('--missmatch', default=2, help='Minimum missmatch in offtargets')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
parser.add_argument('--sg-expected', default=7, help='Expected approximate # of sgRNAs per 100bp. Default: 7')
parser.add_argument('--sg-flex', default=3, help='Room of flexibility for evenness in distribution. For example, 0 denotes perfectly even distribution, 3 means +/- 3bp from perfect evenness is allowed. Default: 3')
parser.add_argument('--plot-off', help="Do not draw plots.", action="store_true")
parser.add_argument('--optimize-off', help="Optimization off.", action="store_true")
parser.add_argument('--distribution-off', help="Distribution filtering off.", action="store_true")
parser.add_argument('--save_tmp', help="Save the temporary files.", action="store_true")
parser.add_argument('--dir', default="tmp", help="Directory to store the figures")
parser.add_argument('--pam-off', action="store_true", help='Do not write PAM sequences to file')

args = parser.parse_args()

#FORMAT = '%(asctime)-15s %(message)s'
FORMAT = ''
log.basicConfig(filename=args.output + '.report.txt', filemode = 'w', format=FORMAT, level=args.verbose)

# set up logging to console
console = log.StreamHandler()
console.setLevel(log.DEBUG)
# add the handler to the root logger
log.getLogger('').addHandler(console)

logger = log.getLogger(__name__)

chrs=[]
for i in range(1,22):
	chrs.append("chr" + str(i))
chrs.extend(["chrX", "chrY"])

sequence_list = []
sequence_list_mutex = threading.Lock()

args.sg_expected = int(args.sg_expected)
args.sg_flex = int(args.sg_flex)

gtf = {}
dhs = {}

sn={'G2': -0.275377128,
	'A3': -0.323887456,
	'C3': 0.172128871,
	'C4': -0.100666209,
	'C5': -0.20180294,
	'G5': 0.245956633,
	'A6': 0.036440041,
	'C6': 0.098376835,
	'C7': -0.741181291,
	'G7': -0.393264397,
	'A12': -0.466099015,
	'A15': 0.085376945,
	'C15': -0.013813972,
	'A16': 0.272620512,
	'C16': -0.119022648,
	'T16': -0.285944222,
	'A17': 0.097454592,
	'G17': -0.17554617,
	'C18': -0.345795451,
	'G18': -0.678096426,
	'A19': 0.22508903,
	'C19': -0.507794051,
	'G20': -0.417373597,
	'T20': -0.054306959,
	'G21': 0.379899366,
	'T21': -0.090712644,
	'C22': 0.057823319,
	'T22': -0.530567296,
	'T23': -0.877007428,
	'C24': -0.876235846,
	'G24': 0.278916259,
	'T24': -0.403102218,
	'A25': -0.077300704,
	'C25': 0.287935617,
	'T25': -0.221637217,
	'G28': -0.689016682,
	'T28': 0.117877577,
	'C29': -0.160445304,
	'G30': 0.386342585
}

dn={'GT2': -0.625778696,
	'GC5': 0.300043317,
	'AA6': -0.834836245,
	'TA6': 0.760627772,
	'GG7': -0.490816749,
	'GG12': -1.516907439,
	'TA12': 0.7092612,
	'TC12': 0.496298609,
	'TT12': -0.586873894,
	'GG13': -0.334563735,
	'GA14': 0.76384993,
	'GC14': -0.53702517,
	'TG17': -0.798146133,
	'GG19': -0.66680873,
	'TC19': 0.353183252,
	'CC20': 0.748072092,
	'TG20': -0.367266772,
	'AC21': 0.568209132,
	'CG21': 0.329072074,
	'GA21': -0.836456755,
	'GG21': -0.782207584,
	'TC22': -1.029692957,
	'CG23': 0.856197823,
	'CT23': -0.463207679,
	'AA24': -0.579492389,
	'AG24': 0.649075537,
	'AG25': -0.077300704,
	'CG25': 0.287935617,
	'TG25': -0.221637217,
	'GT27': 0.117877577,
	'GG29': -0.697740024
}

complement = {'A': 'T',
			  'C': 'G',
			  'G': 'C',
			  'T': 'A',
			  'N': 'N'}

def rev_com(x):
	return ''.join(reversed([complement.get(base,base) for base in list(x)]))

def count_gc(seq):
	seq = seq.upper()
	return seq.count('C') + seq.count('G')

def calculate_gc(seq):
	return count_gc(seq)*100.0/len(seq)

def calculate_efficiency(seq):
	score=0
	t=list(seq)

	gc = count_gc(seq)
	if gc < 10:
		score=0.597636154+(abs(gc-10)*-0.202625894)
	else:
		score=0.597636154+((gc-10)*-0.166587752)

	spacer_length = (30-len(seq))/2 #calculate extended spacer length

	#add single nucleotide score
	for i in range(spacer_length,len(seq)+spacer_length):
		tr=i-spacer_length
		k=t[tr]+str(i)
		score+=sn.get(k,0) #add value from single nucleotide, if not available, add 0
	#add di-nucleotide score
	for i in range(spacer_length,len(seq)-1+spacer_length,2):
		tr=i-spacer_length
		k=t[tr]+t[tr+1]+str(i)
		score+=dn.get(k,0)
	return 1/(1+math.exp(-score))

def write_fasta(dict, dirpath, filename):
	try:
	    os.makedirs(dirpath)
	except OSError:
		pass
	with open(dirpath + '/' + filename, 'w') as f:
		for k in dict:
			f.write('>'+k+'\n'+dict[k]['sequence']+'\n')

def amend_fasta(dict, filename):
	with open(filename, 'a+') as f:
		for k in dict:
			f.write('>'+k+'\n'+dict[k]['sequence']+'\n')

def reset_file(file):
	try:
	    os.remove(file)
	except OSError:
	    pass


def n_lower_chars(string):
    return sum(1 for c in string if c.islower())


def findSg(id, sequence, output_fa):
	sgrnas = {}
	pam=args.PAM.upper()
	l=args.length
	if 'N' in pam:
		l=args.length+pam.count('N') #increase sgRNA length by number of Ns
		pam = pam.replace('N','')

	def findPam(seq,pam,strand):
		sgs = {}
		for idx in range(l,len(seq)-len(pam)):
			if seq[idx:idx+len(pam)].upper() != pam.upper():
				continue
			sgrna_seq = seq[idx-l:idx+len(pam)]
			if 'TTTT' in sgrna_seq.upper():
				continue
			gc = calculate_gc(sgrna_seq)
			if args.gc_min <= gc <= args.gc_max:
				eff = calculate_efficiency(sgrna_seq)
				if eff > 0:
					sgrna_name = id + '_' + str(idx-l) + '_' + strand
					sgs[sgrna_name] = {'pos_start': idx-l,
							   'pos_end': idx + len(pam),
							   'strand': strand,
							   'sequence': sgrna_seq,
							   'efficiency': eff,
							   'GC': gc,
							   'GTF': 0,
							   'DHS': 0,
							   'OTP': 1.0}
		return sgs

	if args.strand == "both" or args.strand == "positive":
		sgrnas = findPam(sequence.upper(), pam, "+")
	if args.strand == "both" or args.strand == "negative":
		sgrnas.update(findPam(rev_com(sequence.upper()), pam, "-"))

	amend_fasta(sgrnas, output_fa)
	return sgrnas


def filter_sg(input_seqs, sgrna_fasta_file):
	cmd = "bowtie -v %d %s %s -f --quiet -p %d"% (args.missmatch - 1, args.bowtie_index, sgrna_fasta_file, args.nthreads)
	bwt_output = subprocess.check_output(shlex.split(cmd))
	bwt_output = bwt_output.strip()
	for match in bwt_output.split('\n'):
		match = match.strip()
		words = match.split('\t')
		ids = words[0].split('_')
		if int(words[6]) > 0:
			try:
				del input_seqs[ids[0]]['sgRNAs'][words[0]]
			except:
				pass


def run_bowtie(fa, output_file):
	cmd = "bowtie %s %s -f -a -v 3 -y --suppress 6,7 --quiet -p %d"% (args.bowtie_index, fa, args.nthreads)
	with open(output_file, 'w') as f:
		subprocess.call(shlex.split(cmd), stdout = f)


def load_tsv(file_path):
	dhs = {}
	with open(file_path) as f:
		for line in f:
			line = line.strip()
			words = line.split('\t')
			if len(words) >= 3:
				dhs.setdefault(words[0],[]).append((int(words[1]), int(words[2])))

	for _,v in dhs.items(): #no use of k, so replaced by underscore
		v.sort(key=lambda tup: tup[0])

	return dhs

def intersect(subject_chr, subject_pos, query_dict):
	first = 0
	last = len(query_dict[subject_chr]) - 1
	while first <= last:
		mid = (first + last)/2
		if subject_pos < query_dict[subject_chr][mid][0]:
			last = mid - 1
		elif subject_pos > query_dict[subject_chr][mid][1]:
			first = mid + 1
		else:
			return True
	return False

def load_bwt(file_path):
	dhs = {}
	with open(file_path) as f:
		for line in f:
			line = line.strip()
			words = line.split('\t')
			if len(words) >= 5:
				dhs.setdefault(words[2],[]).append((int(words[3]), int(words[3])+len(words[4]), words[0]))

	for _,v in dhs.items(): #no use of k, so replaced by underscore
		v.sort(key=lambda tup: tup[0])

	return dhs

def process_bwt(input_seqs, bwt_file, gtf, dhs):
	last_sg = ''
	with open(bwt_file) as f:
		for line in f:
			line = line.strip()
			words = line.split('\t')
			ids = words[0].split('_')
			if words[0] not in input_seqs[ids[0]]['sgRNAs']:
				continue
			if words[2] not in chrs:
				continue
			if words[0] != last_sg:
				perfect_match_counter = 0
			if len(words) < 6:
				perfect_match_counter += 1
				if perfect_match_counter > 1:
					try:
						del input_seqs[ids[0]]['sgRNAs'][words[0]]
					except:
						pass
					continue
			else:
				if words[5].count(':') < args.missmatch:
					try:
						del input_seqs[ids[0]]['sgRNAs'][words[0]]
					except:
						pass
					continue

			if gtf:
				if intersect(str(words[2]),int(words[3]),gtf):
					input_seqs[ids[0]]['sgRNAs'][words[0]]['GTF'] += 1
			if dhs:
				if intersect(str(words[2]),int(words[3]),dhs):
					input_seqs[ids[0]]['sgRNAs'][words[0]]['DHS'] += 1


			last_sg = words[0]

def distribute_sg(input_seqs):

	def update_sg_point(point, k, v):
		sg_point.setdefault(point, {})
		sg_point[point]['sg'] = v
		sg_point[point]['sgid'] = k

	def bresenham(m,n): #Bresenham's line algorithm
		return [i*n//m + n//(2*m) for i in range(m)]

	for seq in input_seqs:
		points = bresenham(args.sg_expected*len(input_seqs[seq]['sequence'])/100, len(input_seqs[seq]['sequence']))
		sg_point = {}
		for k,v in input_seqs[seq]['sgRNAs'].items():
			for point in points:
				if point-args.sg_flex <= v['pos_start'] < point + args.sg_flex:
					if point not in sg_point:
						update_sg_point(point, k, v)
					else:
						if v['OTP'] < sg_point[point]['sg']['OTP']:
							update_sg_point(point, k, v)
						elif v['OTP'] == sg_point[point]['sg']['OTP']:
							if v['efficiency'] > sg_point[point]['sg']['efficiency']:
								update_sg_point(point, k, v)
							elif v['efficiency'] == sg_point[point]['sg']['efficiency']:
								if (abs(point - v['pos_start']) < abs(point - sg_point[point]['sg']['pos_start'])):
									update_sg_point(point, k, v)
					break
		input_seqs[seq]['sgRNAs_filtered'] = {}
		if not args.optimize_off:
			if len(sg_point) < len(points):
				for point in points:
					if point in sg_point:
						continue
					else:
						for k,v in input_seqs[seq]['sgRNAs'].items():
							if point-(100/args.sg_expected) <= v['pos_start'] < point + (100/args.sg_expected):
								update_sg_point(point, k, v)

		for _,v in sg_point.items():
			input_seqs[seq]['sgRNAs_filtered'][v['sgid']] = v['sg']

		pos = {}
		for i in range(len(input_seqs[seq]['sequence'])):
			pos.setdefault(i, 0)
		for _,v in input_seqs[seq]['sgRNAs_filtered'].items():
			for i in range(v['pos_start'] , v['pos_end']):
				pos[i] += 1

		input_seqs[seq]['bp_coverage'] = sum(1 for x in pos.values() if x > 0) / float(len(pos))

		if not args.plot_off:
			x = []
			y = []
			for k,v in pos.items():
				x.append(k)
				y.append(v)
			try:
				os.makedirs(args.dir)
			except OSError:
				pass
			with PdfPages(args.dir + '/' + seq + '.pdf') as pdf:
				plt.figure(figsize=(3, 3))
				plt.plot(x, y)
				plt.title(seq, fontsize = 12)
				plt.xlabel('position', fontsize = 10)
				plt.ylabel('sgRNA coverage', fontsize = 10)
				axes = plt.gca()
				axes.set_ylim([0,max(y) + 1])
				pdf.savefig()
				plt.close()
def test():
	if args.GTF:
		gtf = load_tsv(args.GTF)
		log.info(args.GTF + ' loaded')
	if args.dhs:
		dhs = load_tsv(args.dhs)
		log.info(args.dhs + ' loaded')
		counter=1
		print intersect('chr19', '41197594', gtf)
	with open('/Users/musaahmed/Downloads/tmp2/chr11:2224510-2225010_sgrnas.fa.align.bwt') as f:
		for line in f:
			line = line.strip()
			words = line.split('\t')
			print intersect(str(words[2]), int(words[3]), gtf, dhs)
			if counter == 15:
				break
			counter += 1

def bar_plot(x,y,file_name,ylabel,title):
	with PdfPages(file_name) as pdf:
		plt.figure(figsize=(4, 4))
		plt.bar(x, y, align='center', alpha=0.5)
		plt.ylabel(ylabel)
		plt.xlabel("Input Regions")
		plt.title(title)
		axes = plt.gca()
		axes.set_ylim([0,max(y) * 1.5])
		pdf.savefig()
		plt.close()


def count_total_sg(input_seqs):
	total = 0
	for k in input_seqs:
		total += len(input_seqs[k]['sgRNAs'])
	return total

def count_filtered_sg(input_seqs):
	total = 0
	for k in input_seqs:
		total += len(input_seqs[k]['sgRNAs_filtered'])
	return total


def find_max(dict, feature):
	max = 0
	for _,v in dict['sgRNAs'].items():
		if v[feature] > max:
			max = v[feature]
	return max

def find_otp(input_seqs):
	def calculate_otp(this_gtf, this_dhs, max_gtf, max_dhs):
		if max_gtf > 0:
			otp = int(this_gtf) / float((max_gtf * 1.5))
		else:
			otp = 0
		if max_dhs > 0:
			otp = otp + (int(this_dhs) / float((max_dhs * 3)))
		return otp
	for seq in input_seqs:
 		max_gtf = 0
 		max_dhs = 0
 		if args.GTF:
 			max_gtf = find_max(input_seqs[seq], 'GTF')
		if args.dhs:
			max_dhs = find_max(input_seqs[seq], 'DHS')
 		for k,v in input_seqs[seq]['sgRNAs'].items():
 			input_seqs[seq]['sgRNAs'][k]['OTP'] = calculate_otp(v['GTF'], v['DHS'], max_gtf, max_dhs)

def set_pam_out(v):
	if not args.pam_off:
		write_pam = v['sequence']
	else:
		write_pam = v['sequence'][:-len(args.PAM)]
	return write_pam


def main():
        log.info(args.GTF)
	seq_name = None
	input_seqs = {}
	seq_name = ''
	reset_file(args.output + '_sgrna.fa')
	if args.GTF:
		gtf = load_tsv(args.GTF)
	if args.dhs:
		dhs = load_tsv(args.dhs)
	with open(args.FASTA) as f:
		lines = f.readlines()
		counter = 0
		for line in lines:
			counter += 1
			line = line.strip()
			if line[0]=='>':
				if seq_name:
					input_seqs[seq_name]['sgRNAs']=findSg(seq_name,input_seqs[seq_name]['sequence'],args.output + '_sgrna.fa')
				seq_name = line[1:]
				input_seqs.setdefault(seq_name, {})
				input_seqs[seq_name]['sequence'] = ''
			else:
				input_seqs[seq_name]['sequence'] += line
				if counter == len(lines):
					input_seqs[seq_name]['sgRNAs']=findSg(seq_name,input_seqs[seq_name]['sequence'],args.output + '_sgrna.fa')

	total_sg = count_total_sg(input_seqs)

	log.info("Total sequence: %d"% len(input_seqs))
	log.info("Total sgRNA found: %d"% total_sg)
	filter_sg(input_seqs, args.output +'_sgrna.fa')
	total_sg = count_total_sg(input_seqs)
	log.info("Total sgRNA retained after first pass: %d"% total_sg)
	reset_file(args.output + '_sgrna.filtered.fa')
	for k in input_seqs:
		amend_fasta(input_seqs[k]['sgRNAs'], args.output + '_sgrna.filtered.fa')

	run_bowtie(args.output + '_sgrna.filtered.fa', args.output + '_sgrna.filtered' + '.align.bwt')
	process_bwt(input_seqs, args.output + '_sgrna.filtered' + '.align.bwt', gtf,dhs)
	total_sg = count_total_sg(input_seqs)
	log.info("Total sgRNA retained after second pass: %d"% total_sg)
	find_otp(input_seqs)
	if not args.distribution_off:
		log.info("\nDistribution parameters:\n\tPreferred gap between sgRNAs: %d bp\n\tWiggle room: %d bp\n\tOtimization: %s\n\tPlots: %s" % (100 / args.sg_expected, args.sg_flex, not args.optimize_off, not args.plot_off))
		distribute_sg(input_seqs)
		total_sg = count_filtered_sg(input_seqs)
		log.info("Total sgRNA retained after third pass: %d"% total_sg)
	else:
		log.info("--distribution-off is switched. Distribution optimization skipped.")

	bp_coverage = []
	sg_per_100bp = []
	no_of_sg = []

 	with open(args.output + '.all.txt', 'w') as a, open(args.output + '.sgRNAs.txt', 'w') as f, open(args.output + '.stats.txt', 'w') as l:
	 	l.write("region\ttotal_sg\tretained_sg\tbp_coverage\tno_of_sg\tsg_per_100bp\n")
	 	for seq in input_seqs:
	 		bp_coverage.append(input_seqs[seq]['bp_coverage'])
	 		no_of_sg.append(len(input_seqs[seq]['sgRNAs_filtered']))
	 		sg_per_100bp.append(len(input_seqs[seq]['sgRNAs_filtered']) * 100 / float(len(input_seqs[seq]['sequence'])))
	 		l.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\n"% (seq, len(input_seqs[seq]['sgRNAs']), len(input_seqs[seq]['sgRNAs_filtered']), input_seqs[seq]['bp_coverage'], len(input_seqs[seq]['sgRNAs_filtered']), len(input_seqs[seq]['sgRNAs_filtered']) * 100 / float(len(input_seqs[seq]['sequence']))))
	 		for k,v in input_seqs[seq]['sgRNAs_filtered'].items():
	 			seq_processing = set_pam_out(v)
	 			f.write(seq + '\t' +
	 				str(v['pos_start']) + '\t' +
	 				str(v['pos_end']) + '\t' +
	 				k + '\t' +
	 				seq_processing + '\t' +
	 				"%.2f" % v['efficiency'] + '\t' + "%.2f" % v['OTP'] + '\n')
 			for k,v in input_seqs[seq]['sgRNAs'].items():
 				seq_processing = set_pam_out(v)
	 			a.write(seq + '\t' +
	 				str(v['pos_start']) + '\t' +
	 				str(v['pos_end']) + '\t' +
	 				k + '\t' +
	 				seq_processing + '\t' +
	 				"%.2f" % v['efficiency'] + '\t' + "%.2f" % v['OTP'] + '\n')

	log.info("\nStatistics:\n\tTotal sgRNA: %d\n\tMean sgRNA/region: %.2f\n\tMean sgRNA/100bp: %.2f\n\tRegions with <50%% covered by sgRNAs: %d\n\tAverage coverage (spacer/bp): %.2f" % (count_filtered_sg(input_seqs), sum(no_of_sg)/float(len(no_of_sg)), sum(sg_per_100bp)/float(len(sg_per_100bp)), sum(1 for x in bp_coverage if x < 0.5), sum(bp_coverage)/float(len(bp_coverage))))

	bar_plot(range(1,len(input_seqs)+1), bp_coverage, args.output + '.bp_coverage.pdf', "Coverage", "sgRNA coverage")
	bar_plot(range(1,len(input_seqs)+1), no_of_sg, args.output + '.sgrna_count.pdf', "sgRNA", "No of sgRNAs")

	if not args.save_tmp:
		reset_file(args.output + '_sgrna.fa')
		reset_file(args.output + '_sgrna.filtered.fa')
		reset_file(args.output + '_sgrna.filtered' + '.align.bwt')

if __name__ == "__main__":
	main()
