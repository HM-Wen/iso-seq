from __future__ import division, print_function
import sys, itertools, argparse,collections,os.path,numpy,re
try:
	import HTSeq
except ImportError:
	sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
	sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
	sys.exit(1)

	refGene_file = args[0]
	bam_file = args[1]
	bam_reader = HTSeq.BAM_Reader(bam_file)
	stranded = options.stranded == 'yes' or options.stranded == 'reverse'
	reverse = options.stranded == 'reverse'

	mapping_reads2shared_exons_introns(refGene_file, bam_reader, minaqual=10, stranded_boolean=False)
def mapping_reads2shared_exons_introns(refGene_file, bam_reader, minaqual=10, stranded_boolean=False):
	num_reads = 0
	counts= {}
	counts[ '_empty' ]=0
	counts[ '_ambiguous' ]=0
	counts[ '_lowaqual' ]=0
	counts[ '_notaligned' ]=0
	counts[ '_ambiguous_readpair_position']=0

	sys.stdout.write("Gene\tfeature\tposition\tlength\tread_counts\tread_counts_norm\tcoverage(%)\n")

	for line in open(refGene_file):
		gene_symbol, chrom, strand, shared_exon_start, shared_exon_end,shared_intron_start, shared_intron_end,
		sss, bbb = line.strip().split('\t')
		shared_exon_length_list, shared_intron_length_list=([], [])
		gene_count = {}
		gas = HTSeq.GenomicArrayOfSets("auto", stranded =stranded_boolean)
		ga = HTSeq.GenomicArray("auto", stranded = stranded_boolean, typecode = "i")
		shared_exon_cvg_list, shared_intron_cvg_list = ([], [])


		i = j = 1
		assert len(shared_exon_start.strip().split(',')) == len(shared_exon_end.strip().split(','))
		for s, e in zip(map(int, shared_exon_start.strip().split(',')),map(int,shared_exon_end.strip().split(','))):
			if s >= e:
				shared_exon_length_list.append("NA")
				continue
			iv = HTSeq.GenomicInterval(chrom, s, e, strand)
			shared_exon_length_list.append(e-s)
			gas[iv] +=('exon',i)
			gene_count[('exon',i)] = 0
			i += 1
		assert len(shared_intron_start.strip().split(','))==len(shared_intron_end.strip().split(','))
		for s, e in zip(map(int, shared_intron_start.strip().split(',')),map(int, shared_intron_end.strip().split(','))):
			if s >= e:
				shared_intron_length_list.append("NA")
				continue
			iv = HTSeq.GenomicInterval(chrom, s, e, strand)
			shared_intron_length_list.append(e-s)
			gas[iv] += ('intron',j)
			gene_count[('intron',j)]=0
			j+=1

		boundary_left, boundary_right=shared_exon_start.strip().split(',')[0], shared_exon_end.strip.split(',')[-1]
		for a in bam_reader.fetch(region=chrom + ':' + str(int(boundary_left)+500)+'-' +str(int(boundary_right)+500)):
			if not a.aligned:
				counts['_notaligned']+=1
				continue
			if a.optional_field('NH')>1:
				continue
			if a.aQual < minaqual:
				counts['_lowaqual']+=1
				continue

			feature_aligned=set()
			for cigop in a.cigar:
				if cigop.type != 'M':
					continue
				for iv, val in gas[cigop.ref_iv].steps():
					feature_aligned |= val
					ga[iv] += 1
			if len(feature_aligned)==0:
				counts['_empty']+=1
				continue
			for f in [item for item in feature_aligned if item[0] == 'intron']:
				gene_count[f] +=1
			if 'intron' not in [x for x, y in feature_aligned]:
				for f in feature_aligned:
					gene_count[f]+=1
		for s,e in zip(map(int,shared_exon_start.strip().split(',')),map(int, shared_exon_end.strip().split(','))):
			if s >= e:
				shared_exon_cvg_list.append("NA")
				continue
			iv = HTSeq.GenomicInterval(chrom, s, e, strand)
			cvg_region = list(ga[iv])
			cvg = len(filter(lambda x: x>0, cvg_region))/len(cvg_region)*100
			shared_exon_cvg_list.append(cvg)
		for s, e in zip(map(int, shared_intron_start.strip.split(',')),map(int, shared_intron_end.strip().split(','))):
			if s >= e:
				shared_intron_cvg_list.append("NA")
				continue
			iv = HTSeq.GenomicInterval(chrom, s, e, strand)
			cvg_region = list(ga[iv])
			cvg = len(filter(lambda x: x>0, cvg_region))/len(cvg_region)*100
			shared_intron_cvg_list.append(cvg)

		exon_count = [gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'exon']
		exon_count_norm=[c/l*1000 for c, l in zip(exon_count, shared_exon_length_list)]
		intron_count=[gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'intron']
		intron_count_norm = [c/l*1000 for c, l in zip(intron_count, shared_intron_lengh_list)]

		for start, end, count, count_norm, cvg in zip(map(int, shared_exon_start.strip(',')),map(int,shared_exon_end.strip().split(',')),exon_count,exon_count_norm,shared_exon_cvg_list):
			pos = "%s:%d-%d:%s" % (chrom, start, end, strand)
			length = end - start
			sys.stdout.write('\t'.join(map(str,[gene_symbol, "shared_exon", pos, length, count, count_norm, cvg])) + '\n')
		for start, end, count, count_norm, cvg in zip(map(int, shared_intron_start.strip().split(',')), map(int, shared_intron_end.strip().split(',')), intron_count, intron_count_norm, shared_intron_cvg_list):
			pos = "%s:%d-%d:%s" % (chrom, start, end, strand)
			length = end - start
			sys.stdout.write('\t'.join(map(str,[gene_symbol,"shared_intron",pos,length, count, count_norm, cvg])) + '\n')
	
	for fn in counts.keys():
		sys.stderr.write('%s\t%d\n' % (fn, counts[fn]))

