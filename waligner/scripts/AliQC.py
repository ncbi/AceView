#!/usr/bin/env python3
###############################################################################################
# aliqc.py -- Calculate base alignment statistics for a SAM file
#
# Joe Meehan, NCTR, June 2012
# Danielle et Jean Thierry-Mieg, mieg@ncbi.nlm.nih.gov
#
# Modification history
#
#	8/2012     Convert reads and reference to uppercase, count non-ATGC bases (bad_base_count)
#	9/7/2012   Complement bases on reverse strand for the base matrix counts (A>A,T>A,...)
#	9/8/2012   Count number of alignments per read
#	9/12/2012  Only calculate statistics on the best alignment for a read
#	9/12/2012  Count insert and delete types
#	9/20/2012  Calculate statistics on all alignments, counting proper pairs
#       9/26/2012  mieg: reorder the ATGC in the exported tables of substitutions and indels
#	10/2/2012  Calculate the statistics only for the first alignment of SE, or each end of PE
#       1/10/2012  mieg:Triple delight reorganization of the organigram
#       6/10/2012  mieg:switched tag, now in column 1, and run in the .aliqc.tsv format to simplify typed parsing
#	12/10/2012 jfm: Optimized
###############################################################################################
import getopt
import os
import os.path
import sys
import string
import time
import HTSeq

###############################################################################################
# print usage message and exit
def usage():
	print ("Authors:")
	print ("     Joe Meehan, NCTR, josmeehan@sbcglobal.net")
	print ("     Danielle et Jean Thierry-Mieg, NCBI, mieg@ncbi.mlm.nih.gov")
	print ("     No guarantee whatsoever, no strings attached")
	print ("Usage:")
	print ("     The purpose of this program in phase 1 is to produce QC statistics for SAM/BAM files")
	print ("     then in phase 2 to merge the results, then in phase 3 to format in a user friendly way")
	print ("Mandatory parameters:")
	print ("  --BAM | --SAM | --SAMSORTED | --SAMSORTEDGZ | --SAMGZ | --merge | --view  # select which phase of the program should be run")
	print ("")
	print ("Phase 1:")
	print ("  -S | --SAM       # Phase 1: compute the quality of a SAM file, exports in .aliqc.tsv format")
	print ("  --SAMGZ          # Phase 1: compute the quality of a SAM.gz file, exports in .aliqc.tsv format")
	print ("  --SAMSORTED      # Phase 1: compute the quality of a SAM sorted file, exports in .aliqc.tsv format")
	print ("  --SAMSORTEDGZ    # Phase 1: compute the quality of a SAM gzipped sorted file, exports in .aliqc.tsv format")
	print ("  -B | --BAM       # same for a BAM file, requires `bam2sam` program from Sanger Center")
	print ("  -r r | --run r   # mandatory run identifier")
	print ("  -f f | --fasta f # mandatory the fasta file of the reference used in the BAM/SAM file")
	print ("                   # all results are exported in the 'additive' .aliqc.tsv format")
	print ("   --minAli <int>  # drop SAM lines with less aligned bases, match or mismatch excluding S")
	print ("   -e              # count every base pair, even matches (AA,TT,CC,GG)")
	print ("")
	print ("Phase 2:")
	print ("  -m | --merge     # Phase 2: merge a number of .seqcs.tsv files, output in same format")

	print ("  -r r | --run r   # mandatory run identifier")

# reactivate when the group is working
#	print ("  -r r | --run r   # run identifier, mandatory except if a group file is provided") 
#	print ("  -g g | --group g # The group file 'g' has 2 columns, tab delimited: group_id  run_id")
#	print ("                   # Generate a group count by adding all the run counts")

	print ("                   # The default is to merge all counts into the 'run'")
	print ("                   # all results are exported in the 'additive' .aliqc.tsv format")
	print ("  --nreads int     # optional, number of reads in the fastq file. Defaults to number of reads in the sam file. ")
	print ("  --nbases int     # optional, number of bases in the fastq file. Defaults to number of bases in the sam file.")
	print ("                   # if the BAM file does not include the unaligned reads, a frequent situation,")
	print ("                   # otherwise all percentages are overestimated, giving an artefactual 100/100 aligned")
	print ("")
	print ("Phase 3:")
	print ("  --view html | table     # reads a set of .seqcs.tsv files and export in a user friendly format")
	print ("Options:")
	print ("-i input_file      # optional, input file, defaults to stdin")
	print ("-o file_prefix     # optional, prefix for the output files, the suffix is added by teh program")
	print ("                   # for phase 1 and 2, a single output file is provided, for phase 3, several")
	print ("                   # if -o is not provided, phases 1 and 2 write to stdout")
	print ("                   # in phase 3, the -o argument is mandatory")
	print ("--split            # split the output, rather than exporting a single heterogeneous file")
	print ("")
	print ("Example:")
	print ("     aliqc.py -r Rhs41 -f hs.genome.fasta -B -i f123.bam -o stats/f123  # single run statistics")
	print ("     cat stats/*.aliqc.tsv | aliqc.py --merge -r June -o monthly_stats/june  # merge a set of stats")
	print ("     cat montly_stats/*.aliqc.tsv  | aliqc.py --view table                   # view a global table")
	print ("     cat montly_stats/*.aliqc.tsv  | aliqc.py --view html                    # view as html")
	print ("")
	print ("-? | -h | --help   # This message")
	sys.exit()

###############################################################################################	
# Actual work
# Calculate alignment statistics for a given alignment
#
def calc_seqc( a ):

        # global variables (to accumulate totals across all alignments)
        # global aligned_count # now counted in main loop (jfm 9/28/2012)
        global minAli
        global reads_plus_strand_count
        global reads_minus_strand_count
        global read1_plus_strand_count
        global read1_minus_strand_count
        global read2_plus_strand_count
        global read2_minus_strand_count

        global perfect
        global perfect1
        global perfect2

        global perfect_clipped
        global perfect_clipped1
        global perfect_clipped2

        global intron_N
        global intron_N_1
        global intron_N_2
        
        global intron_n
        global intron_n_1
        global intron_n_2
        
        global prefix_S
        global prefix_S_1
        global prefix_S_2
        
        global suffix_S
        global suffix_S_1
        global suffix_S_2
        
        global prefix_M
        global prefix_M_1
        global prefix_M_2
        
        global suffix_M
        global suffix_M_1
        global suffix_M_2
        
        global prefix_L
        global prefix_L_1
        global prefix_L_2
        
        global suffix_L
        global suffix_L_1
        global suffix_L_2

        global any_substitutions
        global any_substitutions1
        global any_substitutions2

        global any_insertions
        global any_insertions1
        global any_insertions2

        global any_deletions
        global any_deletions1
        global any_deletions2

        global insertions_40
        global insertions_40_1
        global insertions_40_2

        global insertions_1000
        global insertions_1000_1
        global insertions_1000_2

        global insertions_long
        global insertions_long1
        global insertions_long2

        global deletions_40
        global deletions_40_1
        global deletions_40_2

        global deletions_1000
        global deletions_1000_1
        global deletions_1000_2

        global deletions_long
        global deletions_long1
        global deletions_long2

        global sliding_insertions
        global sliding_insertions1
        global sliding_insertions2

        global sliding_deletions
        global sliding_deletions1
        global sliding_deletions2

        global singlet_substitutions
        global singlet_substitutions1
        global singlet_substitutions2

        global singlet_insertions
        global singlet_insertions1
        global singlet_insertions2

        global singlet_deletions
        global singlet_deletions1
        global singlet_deletions2

        global doublet_substitutions
        global doublet_substitutions1
        global doublet_substitutions2

        global doublet_insertions
        global doublet_insertions1
        global doublet_insertions2

        global doublet_deletions
        global doublet_deletions1
        global doublet_deletions2

        global triplet_substitutions
        global triplet_substitutions1
        global triplet_substitutions2

        global longer_substitutions
        global longer_substitutions1
        global longer_substitutions2

        global triplet_insertions
        global triplet_insertions1
        global triplet_insertions2

        global triplet_deletions
        global triplet_deletions1
        global triplet_deletions2

        global quadruplet_substitutions
        global quadruplet_substitutions1
        global quadruplet_substitutions2
        
        global quadruplet_insertions
        global quadruplet_insertions1
        global quadruplet_insertions2
        
        global quadruplet_deletions
        global quadruplet_deletions1
        global quadruplet_deletions2
        
        global quintuplet_substitutions
        global quintuplet_substitutions1
        global quintuplet_substitutions2
        
        global quintuplet_insertions
        global quintuplet_insertions1
        global quintuplet_insertions2
        
        global quintuplet_deletions
        global quintuplet_deletions1
        global quintuplet_deletions2
        
        global sextuplet_substitutions
        global sextuplet_substitutions1
        global sextuplet_substitutions2
        
        global sextuplet_insertions
        global sextuplet_insertions1
        global sextuplet_insertions2
        
        global sextuplet_deletions
        global sextuplet_deletions1
        global sextuplet_deletions2
        
        global septuplet_substitutions
        global septuplet_substitutions1
        global septuplet_substitutions2
        
        global septuplet_insertions
        global septuplet_insertions1
        global septuplet_insertions2
        
        global septuplet_deletions
        global septuplet_deletions1
        global septuplet_deletions2
        
        global octuplet_substitutions
        global octuplet_substitutions1
        global octuplet_substitutions2
        
        global octuplet_insertions
        global octuplet_insertions1
        global octuplet_insertions2
        
        global octuplet_deletions
        global octuplet_deletions1
        global octuplet_deletions2

        global singlet_s_insertions
        global singlet_s_insertions1
        global singlet_s_insertions2

        global singlet_s_deletions
        global singlet_s_deletions1
        global singlet_s_deletions2

        global doublet_s_insertions
        global doublet_s_insertions1
        global doublet_s_insertions2

        global doublet_s_deletions
        global doublet_s_deletions1
        global doublet_s_deletions2

        global triplet_s_insertions
        global triplet_s_insertions1
        global triplet_s_insertions2

        global triplet_s_deletions
        global triplet_s_deletions1
        global triplet_s_deletions2

        global quadruplet_s_insertions
        global quadruplet_s_insertions1
        global quadruplet_s_insertions2
        
        global quadruplet_s_deletions
        global quadruplet_s_deletions1
        global quadruplet_s_deletions2
        
        global quintuplet_s_insertions
        global quintuplet_s_insertions1
        global quintuplet_s_insertions2
        
        global quintuplet_s_deletions
        global quintuplet_s_deletions1
        global quintuplet_s_deletions2
        
        global sextuplet_s_insertions
        global sextuplet_s_insertions1
        global sextuplet_s_insertions2
        
        global sextuplet_s_deletions
        global sextuplet_s_deletions1
        global sextuplet_s_deletions2
        
        global septuplet_s_insertions
        global septuplet_s_insertions1
        global septuplet_s_insertions2
        
        global septuplet_s_deletions
        global septuplet_s_deletions1
        global septuplet_s_deletions2
        
        global octuplet_s_insertions
        global octuplet_s_insertions1
        global octuplet_s_insertions2
        
        global octuplet_s_deletions
        global octuplet_s_deletions1
        global octuplet_s_deletions2
        
        global longer_s_insertions
        global longer_s_insertions1
        global longer_s_insertions2

        global longer_s_deletions
        global longer_s_deletions1
        global longer_s_deletions2

        global match_length_hist
        global match_length_hist1
        global match_length_hist2
        global read_bp
        global read_1_bp
        global read_2_bp
        global ali_bp
        global ali_1_bp
        global ali_2_bp
        global match_base
        global sub_base
        global del_base
        global ins_base
        global match_base_1
        global sub_base_1
        global del_base_1
        global ins_base_1
        global match_base_2
        global sub_base_2
        global del_base_2
        global ins_base_2
        global any_event_hist
        global any_event_hist1
        global any_event_hist2
        global insert_type_hist
        global insert_type_hist1
        global insert_type_hist2
        global delete_type_hist
        global delete_type_hist1
        global delete_type_hist2
        global base_matrix_total
        global base_matrix_total1
        global base_matrix_total2
        global err_pos
        global err_pos1
        global err_pos2
        global err_N1
        global err_N2
        global bad_base_count
        global valid_bases
        global every_base_flag
        global hard_clipped_events
        global hard_clipped_bases
        global soft_clipped_events
        global soft_clipped_bases

        # local variables (for this particular alignment)
        verbose = False
        # jfm 10/5/2012 changed I to uppercase for counting inserts after ref converted to uppercase
        insert_char = "*" 
        soft_char = "*" 
        soft_clip_event = 0
        soft_clip_base = 0
        hard_clip_event = 0
        hard_clip_base = 0
        any_event = 0
        introns = 0
        mapped_ref = ""
        oldInsertPos = -1
        perfect_map = True	# jfm 10/10/2012 Indicates a map with no mismatches, insertions, deletions, introns
        #                Soft Clips are OK.

	# complement and reverse_complement 
        complement = { 'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N' }  # dictionary used to complement bases

        def reverse_complement( seq ):
                rc = ""
                for c in seq[::-1]: 		# step thru seq in reverse order
                        if c in complement :
                                rc += complement[c]	# taking the complement of each base	
                        else:
                                rc += 'N'
                return rc 

        def complement_seq( seq ):
                rc = ""
                for c in seq[::1]: 		# step thru seq in direct order
                        if c in complement :
                                rc += complement[c]	# taking the complement of each base	
                        else:
                                rc += 'N'
                return rc

                
        # reverse complement reads on the reverse strand so they will match the reference
        if a.iv.strand == "-":
                reads_minus_strand_count += 1
                if a.pe_which == "first":
                        read1_minus_strand_count += 1
                elif a.pe_which == "second":
                        read2_minus_strand_count += 1
                mapped_read = str(a.read.get_reverse_complement())  # Use the (possibly optimized) HTSeq reverse complement
        else:
                reads_plus_strand_count += 1
                if a.pe_which == "first":
                        read1_plus_strand_count += 1
                elif a.pe_which == "second":
                        read2_plus_strand_count += 1
                mapped_read = str(a.read)

        ##########
	# Loop through the CIGAR objects to build the mapped reference string (mapped_ref)
        xpos = 0 ;
        for i in range(0, len(a.cigar)):
                type = a.cigar[i].type
                if type == "M" or type == "=" or type == "X": # Alignment match (can be a sequence match or mismatch)
                        mapped_ref += str(seq[a.cigar[i].ref_iv.chrom][ a.cigar[i].ref_iv.start : a.cigar[i].ref_iv.end ])
                elif type == "I": # Insertion to the reference
                        mapped_ref = mapped_ref + ( insert_char * a.cigar[i].size )
                elif type == "S":  # Soft Clip (clipped sequences present in SEQ)                       
                        mapped_ref = mapped_ref + ( soft_char * a.cigar[i].size )
                else : # H/D/n/N/P Deletion or intron or hard-clip or padding
                        continue
		
	# convert mapped_read and mapped_ref to upper case to look for mismatches
        mapped_read = mapped_read.upper()
        mapped_ref = mapped_ref.upper()

        if verbose:
                print ("# mapped_ref*", mapped_ref)
                print ("# mapped_read", mapped_read)
                print ("##", a.read)
                print ()


        ###################
        ##########
	# Loop through the CIGAR objects to get counts 
        for i in range(0, len(a.cigar)):
                type = a.cigar[i].type
                if type == "M" or type == "=" or type == "X": # Alignment match (can be a sequence match or mismatch)
                        xpos += a.cigar[i].size

                elif type == "P": # padding, ignore
                        continue

                elif type == "N": # Intron gt_ag
                        intron_N += 1
                        any_event += 1
                        if a.paired_end:
                                intron_N_1 = 0
                        else:
                                intron_N_2 = 0
                elif type == "n": # Intron gc_ag
                        intron_n += 1
                        any_event += 1
                        if a.paired_end:
                                intron_n_1 = 0
                        else:
                                intron_n_2 = 0

                elif type == "S":
                        # Soft Clip (clipped sequences present in SEQ)
                        iLn = a.cigar[i].size
                        soft_clip_event += 1
                        soft_clip_base += iLn
                        if xpos < 1:
                                if iLn <= 16:
                                        prefix_S += 1
                                elif iLn <= 100:
                                        prefix_M += 1
                                else:
                                        prefix_L += 1
                        else:
                                if iLn <= 16:
                                        suffix_S += 1
                                elif iLn <= 100:
                                        suffix_M += 1
                                else:
                                        suffix_L += 1
                        if a.paired_end:
                                # first and second fragment alignments  
                                if a.pe_which == "first":  
                                        if xpos < 1:
                                                if iLn <= 16:
                                                        prefix_S_1 += 1
                                                elif iLn <= 100:
                                                        prefix_M_1 += 1
                                                else:
                                                        prefix_L_1 += 1
                                        else:
                                                if iLn <= 16:
                                                        suffix_S_1 += 1
                                                elif iLn <= 100:
                                                        suffix_M_1 += 1
                                                else:
                                                        suffix_L_1 += 1
                                else:
                                        if xpos < 1:
                                                if iLn <= 16:
                                                        prefix_S_2 += 1
                                                elif iLn <= 100:
                                                        prefix_M_2 += 1
                                                else:
                                                        prefix_L_2 += 1
                                        else:
                                                if iLn <= 16:
                                                        suffix_S_2 += 1
                                                elif iLn <= 100:
                                                        suffix_M_2 += 1
                                                else:
                                                        suffix_L_2 += 1
                        xpos += iLn


                elif type == "H": 
                        # hard Clip (clipped sequences absent in seq SEQ)
                        iLn = a.cigar[i].size
                        hard_clip_event += 1
                        hard_clip_base += iLn
                        


                elif type == "D": # Deletion from the reference
                        perfect_map = False # jfm 10/10/2012
                        iLn = a.cigar[i].size

                        any_event += 1
                        any_deletions += 1 
                        if iLn <= 9:
                                del_base += iLn
                        # count deletion types
                        deleton = str(seq[a.cigar[i].ref_iv.chrom][ a.cigar[i].ref_iv.start : a.cigar[i].ref_iv.end])
                        deleton = deleton.upper()  
                        if iLn < 9:
                                if a.iv.strand == "-": 	
                                        # complement deleted bases on the minus strand
                                        k = reverse_complement(deleton)
                                        # to reflect what should have been read
                                else:
                                        k = deleton
                                if k in delete_type_hist :
                                                delete_type_hist[k] += 1
                                else:
                                        delete_type_hist[k] = 1
                                if a.paired_end:
                                        # first and second fragment alignments  
                                        if a.pe_which == "first":
                                                del_base_1 += iLn
                                                if k in delete_type_hist1 :
                                                        delete_type_hist1[k] += 1
                                                else:
                                                        delete_type_hist1[k] = 1
                                        elif a.pe_which == "second":
                                                del_base_2 += iLn
                                                if k in delete_type_hist2 :
                                                        delete_type_hist2[k] += 1
                                                else:
                                                        delete_type_hist2[k] = 1
 
                                # record position of very short insertions (sequencing error spikes)
                                pos = xpos
                                if a.iv.strand == "-":  # unless we're on the minus strand 
                                        pos = len(mapped_read) - pos - 1 - 1 # base 0 of 100 becomes (100 - 0 - 1) = 99, zero based coordinates
                                if pos in err_pos :
                                        err_pos[pos] += 1  # mieg: single event, do not count len(delstr)
                                else:
                                        err_pos[pos] = 1
                                if a.paired_end:
                                        # first and second fragment alignments  
                                        if a.pe_which == "first":
                                                # Count mismatches at this position
                                                if pos in err_pos1:
                                                        err_pos1[pos] += 1
                                                else:
                                                        err_pos1[pos] = 1
                                        elif a.pe_which == "second":
                                                # Count mismatches at this position
                                                if pos in err_pos2:
                                                        err_pos2[pos] += 1
                                                else:
                                                        err_pos2[pos] = 1



                        # evaluate sliding
                        if a.cigar[i].ref_iv.start < iLn or a.cigar[i].ref_iv.start + 2 * iLn > len(seq[a.cigar[i].ref_iv.chrom]):
                                iLn = 0
                        deleton = str(seq[a.cigar[i].ref_iv.chrom][ a.cigar[i].ref_iv.start - iLn : a.cigar[i].ref_iv.end + iLn ])
                        deleton = deleton.upper()    
                        # Count singlet, doublet, triplet with or without sliding
                        isSliding = 0
                        iLnS = 0
                        for j in range (iLn):
                                if deleton[iLn + j] == deleton[2*iLn+j]:
                                        iLnS += 1
                                else:
                                        break
                        for j in range (iLn):
                                if deleton[2*iLn -1 - j] == deleton[iLn - 1 - j]:
                                        iLnS += 1
                                else:
                                        break
                        if iLnS >= iLn:
                                if iLn <= 8:
                                        iLnS = - iLn   # sliding
                                sliding_insertions += 1 
                                isSliding = 1
                        # print (deleton,iLn,iLnS, isSliding)
                        if iLn == 1:
                                singlet_deletions += 1
                        elif iLn == 2:
                                doublet_deletions += 1
                        elif iLn == 3:
                                triplet_deletions += 1
                        elif iLn == 4:
                                quadruplet_deletions += 1
                        elif iLn == 5:
                                quintuplet_deletions += 1
                        elif iLn == 6:
                                sextuplet_deletions += 1
                        elif iLn == 7:
                                septuplet_deletions += 1
                        elif iLn == 8:
                                octuplet_deletions += 1
                        elif iLn < 41:
                                deletions_40 += 1
                        elif iLn < 1001:
                                deletions_1000 += 1
                        else:
                                deletions_long += 1
                                

                        if iLnS == -1:
                                singlet_s_deletions += 1
                        elif iLnS == -2:
                                doublet_s_deletions += 1
                        elif iLnS == -3:
                                triplet_s_deletions += 1
                        elif iLnS == -4:
                                quadruplet_s_deletions += 1
                        elif iLnS == -5:
                                quintuplet_s_deletions += 1
                        elif iLnS == -6:
                                sextuplet_s_deletions += 1
                        elif iLnS == -7:
                                septuplet_s_deletions += 1
                        elif iLnS == -8:
                                octuplet_s_deletions += 1

                        if a.paired_end:
                                # first and second fragment alignments  
                                if a.pe_which == "first":  
                                        any_deletions1 += 1
                                        if isSliding:
                                                sliding_deletions1 += 1
                                                if iLn == 1:
                                                        singlet_deletions1 += 1
                                                elif iLn == 2:
                                                        doublet_deletions1 += 1
                                                elif iLn == 3:
                                                        triplet_deletions1 += 1
                                                elif iLn == 4:
                                                        quadruplet_deletions1 += 1
                                                elif iLn == 5:
                                                        quintuplet_deletions1 += 1
                                                elif iLn == 6:
                                                        sextuplet_deletions1 += 1
                                                elif iLn == 7:
                                                        septuplet_deletions1 += 1
                                                elif iLn == 8:
                                                        octuplet_deletions1 += 1
                                                elif iLn < 41:
                                                        deletions_40_1 += 1
                                                elif iLn < 1001:
                                                        deletions_1000_1 += 1
                                                else:
                                                        deletions_long1 += 1
                                                        
                                                if iLnS == -1:
                                                        singlet_s_deletions1 += 1
                                                elif iLnS == -2:
                                                        doublet_s_deletions1 += 1
                                                elif iLnS == -3:
                                                        triplet_s_deletions1 += 1
                                                elif iLnS == -4:
                                                        quadruplet_s_deletions1 += 1
                                                elif iLnS == -5:
                                                        quintuplet_s_deletions1 += 1
                                                elif iLnS == -6:
                                                        sextuplet_s_deletions1 += 1
                                                elif iLnS == -7:
                                                        septuplet_s_deletions1 += 1
                                                elif iLnS == -8:
                                                        octuplet_s_deletions += 1

                                        else:
                                                any_deletions2 += 1
                                                if isSliding:
                                                        sliding_deletions2 += 1
                                                if iLn == 1:
                                                        singlet_deletions2 += 1
                                                elif iLn == 2:
                                                        doublet_deletions2 += 1
                                                elif iLn == 3:
                                                        triplet_deletions2 += 1
                                                elif iLn == 4:
                                                        quadruplet_deletions2 += 1
                                                elif iLn == 5:
                                                        quintuplet_deletions2 += 1
                                                elif iLn == 6:
                                                        sextuplet_deletions2 += 1
                                                elif iLn == 7:
                                                        septuplet_deletions2 += 1
                                                elif iLn == 8:
                                                        octuplet_deletions2 += 1
                                                elif iLn < 41:
                                                        deletions_40_2 += 1
                                                elif iLn < 1001:
                                                        deletions_1000_2 += 1
                                                else:
                                                        deletions_long2 += 1
                                                        
                                                if iLnS == -1:
                                                        singlet_s_deletions2 += 1
                                                elif iLnS == -2:
                                                        doublet_s_deletions2 += 1
                                                elif iLnS == -3:
                                                        triplet_s_deletions2 += 1
                                                elif iLnS == -4:
                                                        quadruplet_s_deletions2 += 1
                                                elif iLnS == -5:
                                                        quintuplet_s_deletions2 += 1
                                                elif iLnS == -6:
                                                        sextuplet_s_deletions2 += 1
                                                elif iLnS == -7:
                                                        septuplet_s_deletions2 += 1
                                                elif iLnS == -8:
                                                        octuplet_s_deletions2 += 1
                                
                elif type == "I": # Insertion to the reference
                        perfect_map = False 
                        iLn = a.cigar[i].size 
                        xpos2 = xpos + iLn

                        any_event += 1
                        any_insertions += 1 
                        if iLn <= 8:
                                ins_base += iLn
                        # count insertion types
                        inseron = str(mapped_read)[xpos : xpos2 ]
                        inseron = inseron.upper()   
                        # print (xpos, inseron,str(mapped_read)[xpos: xpos+10 ])
                        if iLn <= 8:
                                if a.iv.strand == "-": 	
                                        # complement deleted bases on the minus strand
                                        k = reverse_complement(inseron)
                                        # to reflect what should have been read
                                else:
                                        k = inseron
                                if k in insert_type_hist :
                                                insert_type_hist[k] += 1
                                else:
                                        insert_type_hist[k] = 1
                                if a.paired_end:
                                        # first and second fragment alignments  
                                        if a.pe_which == "first":
                                                ins_base_1 += iLn
                                                if k in insert_type_hist1 :
                                                        insert_type_hist1[k] += 1
                                                else:
                                                        insert_type_hist1[k] = 1
                                        elif a.pe_which == "second":
                                                ins_base_2 += iLn
                                                if k in insert_type_hist2 :
                                                        insert_type_hist2[k] += 1
                                                else:
                                                        insert_type_hist2[k] = 1
 
                                # record position of very short insertions (sequencing error spikes)
                                pos = xpos
                                if a.iv.strand == "-":  # unless we're on the minus strand 
                                        pos = len(mapped_read) - pos - 1 - 1 # base 0 of 100 becomes (100 - 0 - 1) = 99, zero based coordinates
                                if pos in err_pos :
                                        err_pos[pos] += 1  # mieg: single event, do not count len(delstr)
                                else:
                                        err_pos[pos] = 1
                                if a.paired_end:
                                        # first and second fragment alignments  
                                        if a.pe_which == "first":
                                                # Count mismatches at this position
                                                if pos in err_pos1:
                                                        err_pos1[pos] += 1
                                                else:
                                                        err_pos1[pos] = 1
                                        elif a.pe_which == "second":
                                                # Count mismatches at this position
                                                if pos in err_pos2:
                                                        err_pos2[pos] += 1
                                                else:
                                                        err_pos2[pos] = 1

                        # evaluate sliding, 
                        if xpos < iLn or xpos + 2 * iLn > len(mapped_read):
                                iLn = 0
                        inseron = str(mapped_read)[xpos - iLn : xpos2 + iLn ]
                        inseron = inseron.upper()    
                        # Count singlet, doublet, triplet with or without sliding
                        isSliding = 0
                        iLnS = 0
                        for j in range (iLn):
                                if inseron[iLn + j] == inseron[2*iLn+j]:
                                        iLnS += 1
                                else:
                                        break
                        for j in range (iLn):
                                if inseron[2*iLn -1 - j] == inseron[iLn - 1 - j]:
                                        iLnS += 1
                                else:
                                        break
                        if iLnS >= iLn:
                                if iLn <= 8:
                                        iLnS = - iLn   # sliding
                                sliding_insertions += 1 
                                isSliding = 1
                        if iLn == 1:
                                singlet_insertions += 1
                        elif iLn == 2:
                                doublet_insertions += 1
                        elif iLn == 3:
                                triplet_insertions += 1
                        elif iLn == 4:
                                quadruplet_insertions += 1
                        elif iLn == 5:
                                quintuplet_insertions += 1
                        elif iLn == 6:
                                sextuplet_insertions += 1
                        elif iLn == 7:
                                septuplet_insertions += 1
                        elif iLn == 8:
                                octuplet_insertions += 1
                        elif iLn < 41:
                                insertions_40 += 1
                        elif iLn < 1001:
                                insertions_1000 += 1
                        else:
                                insertions_long += 1
                                

                        if iLnS == -1:
                                singlet_s_insertions += 1
                        elif iLnS == -2:
                                doublet_s_insertions += 1
                        elif iLnS == -3:
                                triplet_s_insertions += 1
                        elif iLnS == -4:
                                quadruplet_s_insertions += 1
                        elif iLnS == -5:
                                quintuplet_s_insertions += 1
                        elif iLnS == -6:
                                sextuplet_s_insertions += 1
                        elif iLnS == -7:
                                septuplet_s_insertions += 1
                        elif iLnS == -8:
                                octuplet_s_insertions += 1

                        if a.paired_end:
                                # first and second fragment alignments  
                                if a.pe_which == "first":  
                                        any_insertions1 += 1
                                        if isSliding:
                                                sliding_insertions1 += 1
                                                if iLn == 1:
                                                        singlet_insertions1 += 1
                                                elif iLn == 2:
                                                        doublet_insertions1 += 1
                                                elif iLn == 3:
                                                        triplet_insertions1 += 1
                                                elif iLn == 4:
                                                        quadruplet_insertions1 += 1
                                                elif iLn == 5:
                                                        quintuplet_insertions1 += 1
                                                elif iLn == 6:
                                                        sextuplet_insertions1 += 1
                                                elif iLn == 7:
                                                        septuplet_insertions1 += 1
                                                elif iLn == 8:
                                                        octuplet_insertions1 += 1
                                                elif iLn < 41:
                                                        insertions_40_1 += 1
                                                elif iLn < 1001:
                                                        insertions_1000_1 += 1
                                                else:
                                                        insertions_long1 += 1
                                                        
                                                if iLnS == -1:
                                                        singlet_s_insertions1 += 1
                                                elif iLnS == -2:
                                                        doublet_s_insertions1 += 1
                                                elif iLnS == -3:
                                                        triplet_s_insertions1 += 1
                                                elif iLnS == -4:
                                                        quadruplet_s_insertions1 += 1
                                                elif iLnS == -5:
                                                        quintuplet_s_insertions1 += 1
                                                elif iLnS == -6:
                                                        sextuplet_s_insertions1 += 1
                                                elif iLnS == -7:
                                                        septuplet_s_insertions1 += 1
                                                elif iLnS == -8:
                                                        octuplet_s_insertions += 1

                                        else:
                                                any_insertions2 += 1
                                                if isSliding:
                                                        sliding_insertions2 += 1
                                                if iLn == 1:
                                                        singlet_insertions2 += 1
                                                elif iLn == 2:
                                                        doublet_insertions2 += 1
                                                elif iLn == 3:
                                                        triplet_insertions2 += 1
                                                elif iLn == 4:
                                                        quadruplet_insertions2 += 1
                                                elif iLn == 5:
                                                        quintuplet_insertions2 += 1
                                                elif iLn == 6:
                                                        sextuplet_insertions2 += 1
                                                elif iLn == 7:
                                                        septuplet_insertions2 += 1
                                                elif iLn == 8:
                                                        octuplet_insertions2 += 1
                                                elif iLn < 41:
                                                        insertions_40_2 += 1
                                                elif iLn < 1001:
                                                        insertions_1000_2 += 1
                                                else:
                                                        insertions_long2 += 1
                                                        
                                                if iLnS == -1:
                                                        singlet_s_insertions2 += 1
                                                elif iLnS == -2:
                                                        doublet_s_insertions2 += 1
                                                elif iLnS == -3:
                                                        triplet_s_insertions2 += 1
                                                elif iLnS == -4:
                                                        quadruplet_s_insertions2 += 1
                                                elif iLnS == -5:
                                                        quintuplet_s_insertions2 += 1
                                                elif iLnS == -6:
                                                        sextuplet_s_insertions2 += 1
                                                elif iLnS == -7:
                                                        septuplet_s_insertions2 += 1
                                                elif iLnS == -8:
                                                        octuplet_s_insertions2 += 1

                        xpos += a.cigar[i].size

                else:
                        # Something completely different
                        perfect_map = False 
                        print ("# Unexpected type:",a.cigar[i].type, a.cigar[i].size, a.cigar[i].ref_iv.chrom,a.cigar[i].ref_iv.strand)


        # Update aligned length histogram (Total of CIGAR M sizes)

        # k = aligned_length
        k = len(mapped_read) - soft_clip_base
        k1 = len(mapped_read) + hard_clip_base
        if (k < minAli):
                return
        # print ("Len", k)
        ali_bp += k
        read_bp += k1

        if a.paired_end:
                # first and second fragment alignments  
                if a.pe_which == "first":
                        ali_1_bp +=  k
                        read_1_bp += k1
                elif a.pe_which == "second":
                        ali_2_bp +=  k
                        read_2_bp += k1

        ###################

        # Characterize substitutions, but protect them by 4 exact bases left and right since long X ranges are mushy
        if  len(mapped_read) != len(mapped_ref) :
                print ("mapped_read=",mapped_read,"  ln=", len(mapped_read)) ;
                print ("mapped_ref =",mapped_ref, "   ln=",len(mapped_ref)) ;
                print (a)
                return
        assert len(mapped_read) == len(mapped_ref)
        m_base = 0
        if mapped_read == mapped_ref:
                m_base = k
        else:
                kk = ["a","a","a","a","a","a","a","a","a","a","a","a","a"]
                # print (kk)
                inSub = 0
                presub = 0
                postsub1 = 0
                postsub2 = 0
                postsub3 = 0
                postsub4 = 0
                postsub5 = 0
                postsub6 = 0
                postsub7 = 0
                postsub8 = 0
                iMax = len(mapped_read)
                for i in range( len(mapped_read) ):

                        # Determine position of base in physical read (jfm 10/4/2012)
                        pos = i
                        if a.iv.strand == "-":
                                pos = len(mapped_read) - i - 1 # base 0 of 100 becomes (100 - 0 - 1) = 99, zero based coordinates

                        # Count N's at this position (jfm 10/4/2012) 
                        if mapped_read[i] == 'N':
                                if a.paired_end:
                                        if a.pe_which == "first":
                                                if pos in err_N1:
                                                        err_N1[pos] += 1
                                                else:
                                                        err_N1[pos] = 1
                                        elif a.pe_which == "second":
                                                if pos in err_N2 :
                                                        err_N2[pos] += 1
                                                else:
                                                        err_N2[pos] = 1

                        if (mapped_read[i] in valid_bases) and (mapped_ref[i] in valid_bases):

                                # Count mismatches and error positions in this read
                                if i >= iMax - 2 or mapped_read[i+2] != mapped_ref[i+2]:
                                        presub = 0
                                else:
                                        presub += 1

                                if i < 2 or mapped_read[i-2] != mapped_ref[i-2]:
                                        postsub1 = 0
                                else:
                                        postsub1 += 1
                                if i < 3 or mapped_read[i-3] != mapped_ref[i-3]:
                                        postsub2 = 0
                                else:
                                        postsub2 += 1
                                if i < 4 or mapped_read[i-4] != mapped_ref[i-4]:
                                        postsub3 = 0
                                else:
                                        postsub3 += 1
                                if i < 5 or mapped_read[i-5] != mapped_ref[i-5]:
                                        postsub4 = 0
                                else:
                                        postsub4 += 1
                                if i < 6 or mapped_read[i-6] != mapped_ref[i-6]:
                                        postsub5 = 0
                                else:
                                        postsub5 += 1
                                if i < 7 or mapped_read[i-7] != mapped_ref[i-7]:
                                        postsub6 = 0
                                else:
                                        postsub6 += 1
                                if i < 8 or mapped_read[i-8] != mapped_ref[i-8]:
                                        postsub7 = 0
                                else:
                                        postsub7 += 1
                                if i < 9 or mapped_read[i-9] != mapped_ref[i-9]:
                                        postsub8 = 0
                                else:
                                        postsub8 += 1

                                if (mapped_read[i] != mapped_ref[i]):
                                        if (inSub < 8):
                                                # Build base matrix key = base in reference + base in read
                                                # reverse complement the bases if they are on the reverse strand (9/7/2012)
                                                if a.iv.strand == "-":
                                                        kk[inSub] = complement[mapped_ref[i]] + complement[mapped_read[i]]
                                                else:
                                                        kk[inSub] = mapped_ref[i] + mapped_read[i]
                                                #print (inSub, kk, presub, postsub1, postsub2, postsub3)
                                        inSub += 1

                                else:
                                        m_base += 1
                                        if inSub > 0:
                                                any_substitutions += 1
                                                sub_base += inSub
                                                perfect_map = False # jfm 10/10/2012
                                                any_event += 1
                                                # print (inSub, kk, presub, postsub1, postsub2, postsub3)
                                        if inSub == 1:
                                                singlet_substitutions += 1
                                        elif inSub == 2:
                                                doublet_substitutions += 1
                                        elif inSub == 3:
                                                triplet_substitutions += 1
                                        elif inSub == 4:
                                                quadruplet_substitutions += 1
                                        elif inSub == 5:
                                                quintuplet_substitutions += 1
                                        elif inSub == 6:
                                                sextuplet_substitutions += 1
                                        elif inSub == 7:
                                                septuplet_substitutions += 1
                                        elif inSub == 8:
                                                octuplet_substitutions += 1
                                        elif inSub > 8:
                                                longer_substitutions += 1
                                        if a.paired_end:
                                                # first and second fragment alignments  
                                                if a.pe_which == "first":
                                                        if inSub > 0:
                                                                sub_base_1 += inSub
                                                                any_substitutions1 += 1
                                                                if inSub == 1:
                                                                        singlet_substitutions1 += 1
                                                                elif inSub == 2:
                                                                        doublet_substitutions1 += 1
                                                                elif inSub == 3:
                                                                        triplet_substitutions1 += 1
                                                                elif inSub == 4:
                                                                        quadruplet_substitutions1 += 1
                                                                elif inSub == 5:
                                                                        quintuplet_substitutions1 += 1
                                                                elif inSub == 6:
                                                                        sextuplet_substitutions1 += 1
                                                                elif inSub == 7:
                                                                        septuplet_substitutions1 += 1
                                                                elif inSub == 8:
                                                                        octuplet_substitutions1 += 1
                                                                elif inSub > 8:
                                                                        longer_substitutions1 += 1
 
                                                elif a.pe_which == "second":
                                                        if inSub > 0:
                                                                sub_base_2 += inSub
                                                                any_substitutions2 += 1
                                                                if inSub == 1:
                                                                        singlet_substitutions2 += 1
                                                                elif inSub == 2:
                                                                        doublet_substitutions2 += 1
                                                                elif inSub == 3:
                                                                        triplet_substitutions2 += 1
                                                                elif inSub == 4:
                                                                        quadruplet_substitutions2 += 1
                                                                elif inSub == 5:
                                                                        quintuplet_substitutions2 += 1
                                                                elif inSub == 6:
                                                                        sextuplet_substitutions2 += 1
                                                                elif inSub == 7:
                                                                        septuplet_substitutions2 += 1
                                                                elif inSub == 8:
                                                                        octuplet_substitutions2 += 1
                                                                elif inSub > 8:
                                                                        longer_substitutions2 += 1

                                        # report the nature of the substitution only in non fuzzy zones
                                        if inSub == 1 and (presub < 3 or postsub1 < 3):
                                                kk[0] = 0
                                        if inSub == 2 and (presub < 3 or postsub2 < 3):
                                                kk[1] = 0
                                        if inSub == 3 and (presub < 3 or postsub3 < 3):
                                                kk[2] = 0
                                        if inSub == 4 and (presub < 3 or postsub4 < 3):
                                                kk[2] = 0
                                        if inSub == 5 and (presub < 3 or postsub5 < 3):
                                                kk[2] = 0
                                        if inSub == 6 and (presub < 3 or postsub6 < 3):
                                                kk[2] = 0
                                        if inSub == 7 and (presub < 3 or postsub7 < 3):
                                                kk[2] = 0
                                        if inSub == 8 and (presub < 3 or postsub8 < 3):
                                                kk[2] = 0
                                        if inSub <= 8:
                                                for ik in range (inSub):
                                                        k = kk[ik]
                                                        if a.iv.strand == "-":
                                                                posx = pos + 1 + ik
                                                        else:
                                                                posx = pos - 1 - ik
                                                        #  print ("---",pos, ik, k, presub,postsub3)
                                                        if k == 0:
                                                                continue
                                                        # count substitution types
                                                        base_matrix_total[k] += 1
                                                        # Count mismatches at this position
                                                        if posx in err_pos:
                                                                err_pos[posx] += 1
                                                        else:
                                                                err_pos[posx] = 1
                                                
                                                        # first and second fragment alignments  
                                                        if a.paired_end:
                                                                if a.pe_which == "first":
                                                                        base_matrix_total1[k] += 1
                                                                        if posx in err_pos1:
                                                                                err_pos1[posx] += 1
                                                                        else:
                                                                                err_pos1[posx] = 1
                                                                elif a.pe_which == "second":
                                                                        base_matrix_total2[k] += 1
                                                                        if posx in err_pos2:
                                                                                err_pos2[posx] += 1
                                                                        else:
                                                                                err_pos2[posx] = 1


                                        inSub = 0
                        else:
                                bad_base_count += 1

        match_base += m_base
        k = m_base
        if k in match_length_hist :
                match_length_hist[k] += 1
        else:
                match_length_hist[k] = 1
        if a.paired_end:
                # first and second fragment alignments  
                if a.pe_which == "first":
                        match_base_1 += m_base 
                        if k in match_length_hist1 :
                                match_length_hist1[k] += 1
                        else:
                                match_length_hist1[k] = 1
                elif a.pe_which == "second":
                        match_base_2 +=  k
                        if k in match_length_hist2 :
                                match_length_hist2[k] += 1
                        else:
                                match_length_hist2[k] = 1

        # Update mismatch rate histogram and perfect map counters
        soft_clipped_bases += soft_clip_base
        soft_clipped_events += soft_clip_event
        hard_clipped_bases += hard_clip_base
        hard_clipped_events += hard_clip_event
        if perfect_map:
                if soft_clip_base == 0:
                        perfect += 1
                else:
                        perfect_clipped += 1
        k = any_event
        if k in any_event_hist:
                any_event_hist[k] += 1
        else:
                any_event_hist[k] = 1
        if a.paired_end:
                # first and second fragment alignments  
                if a.pe_which == "first":
                        if perfect_map:
                                if soft_clip_event == 0:
                                        perfect1 += 1
                                else:
                                        perfect_clipped1 += 1
                        if k in any_event_hist1:
                                any_event_hist1[k] += 1
                        else:
                                any_event_hist1[k] = 1
                elif a.pe_which == "second":
                        if perfect_map:
                                if soft_clip_event == 0:
                                        perfect2 += 1
                                else:
                                        perfect_clipped2 += 1
                        if k in any_event_hist2:
                                any_event_hist2[k] += 1
                        else:
                                any_event_hist2[k] = 1

        return 1

###############################################################################################	
###############################################################################################	
# Phase2 merge the group counts and reexport
# Phase2a: parse the input, recognize the runId in column 1. 

def phase2Parse( f,run,tags,runs,tag_runs,tags2,tag_val,various,merge) :
        n = 0
        various["minAliLn"]=0
        various["maxAliLn"]=0
        various["maxError"]=0
        various["maxErrPos"]=0
        various["maxErrPos1"]=0
        various["maxErrPos2"]=0
        various["maxMultiAli"]=1
        words = {}
        try:
                for line in f:
                        n = n + 1
                        words=line.split("\t")
                        tag=words[0]
                        if len(words) < 2 or tag[0] == "#":
                                continue
                        if tag[0:9] == "MultiAli:" :
                                x = int (tag[9:len(tag)])
                                if (x > various["maxMultiAli"]):
                                        various["maxMultiAli"] = x
                        if tag[0:16] == "MismatchPerRead:" :
                                x = int (tag[16:len(tag)])
                                if (x > 100):
                                        x = 100
                                if (x > various["maxError"]):
                                        various["maxError"] = x	
                        if tag[0:14] == "Ambiguous_pos:" :
                                x = int (tag[14:len(tag)])
                                if (x > 200):
                                        x = 200
                                tag = "Ambiguous_pos:%04d" % (x)
                        if tag[0:8] == "Err_pos:" :
                                x = int (tag[8:len(tag)])
                                if (x > 200):
                                        x = 200
                                tag = "Err_pos:%04d" % (x)
                                if (len(words) > 3 and int(words[3]) > 0 and x > various["maxErrPos"]):
                                        various["maxErrPos"] = x
                                if (len(words) > 4 and int(words[4]) > 0 and x > various["maxErrPos1"]):
                                        various["maxErrPos1"] = x
                                if (len(words) > 5 and int(words[5]) > 0 and x > various["maxErrPos2"]):
                                        various["maxErrPos2"] = x
                        if tag[0:6] == "AliLn:" :
                                x = int (tag[6:len(tag)])
                                if (x > 200):
                                        x = 200
                                tag = "AliLn:%04d" % (x)
                                if (x > various["maxAliLn"]):
                                        various["maxAliLn"] = x
                                if (x < various["minAliLn"] or  various["minAliLn"] == 0):
                                        various["minAliLn"] = x
                        if merge == 0 :
                                # ignore the given name, we ant to merge
                                run=words[1]
                        # print ("...",words[0],"---",words[1],"===",words[2])
                        nn=int(words[2])
                        tag_run = tag+tab+run
                        tags[tag] = 1
                        runs[run] = 1
                        tag_runs[tag_run] = nn
                        j = 0
                        while (j < nn and j + 3 < len(words)) :
                                j = j + 1
                                x = int(float(words[j + 2]))
                                tag2 = tag_run + "#" + str(j)
                                if tag2 in tag_val:
                                        tag_val[tag2] += x 
                                else:
                                        tag_val[tag2] = x 

        finally:
                f.close()

###############################################################################################	
# reexport the cumulated counts

def phase2Export (f,runs, tags,tag_runs,tags2,tag_val,nReadsGiven,nBasesGiven ) :
        for run in sorted(runs) :
                tag2 = "Reads_in_file" + tab + run + "#3"
                nrf = 0
                if tag2 in tag_val:
                        nrf = tag_val[tag2]
                        
                for tag in tags:
                        if nReadsGiven and tag == "Reads_in_run" :
                                continue
                        tag_run = tag + tab + run
                        if tag_run in tag_runs :
                                nn = tag_runs[tag_run]
                                # print ("tag_run=",tag_run," nn=",nn)
                                j = 0
                                f.write (tag_run+tab +str(nn))
                                if nrf == 0 and tag == "Reads_in_run" :
                                        tag2 = tag_run + "#1" ;
                                        if tag2 in tag_val:
                                                x = tag_val[tag2]
                                                tag_val[ tag_run + "#2"] = x
                                                tag_val[ tag_run + "#3"] = 0
                                while(j < nn) :
                                        j = j + 1
                                        tag2 = tag_run + "#" + str(j)
                                        # print ("nn=",nn," j=",j," tag2=",tag2)
                                        if tag2 in tag_val:
                                                f.write (tab+str(tag_val[tag2]))
                        f.write(newline)
                if nReadsGiven > 0 :
                        if nrf > 0 :
                                f.write ("Reads_in_run"+tab+run_id+"\t3"+tab+str(nReadsGiven)+tab+str(nReadsGiven/2)+tab+str(nReadsGiven/2)+newline)
                        else:
                                f.write ("Reads_in_run"+tab+run_id+"\t3"+tab+str(nReadsGiven)+tab+str(nReadsGiven)+tab+str(0)+newline)
                if nBasesGiven > 0 :
                        if nrf > 0 :
                                f.write ("Bases_in_run"+tab+run_id+"\t3"+tab+str(nBasesGiven)+tab+str(nBasesGiven/2)+tab+str(nBasesGiven/2)+newline)
                        else:
                                f.write ("Bases_in_run"+tab+run_id+"\t3"+tab+str(nBasesGiven)+tab+str(nBasesGiven)+tab+str(0)+newline)


###############################################################################################	
# Phase2 merge the group counts and reexport

def phase2(run, input_stream, tsv_file,nReadsGiven,nBasesGiven  ):
	tags = {}
	tags2 = {}
	runs = {}
	tag_val = {}
	tag_runs = {}
	tag_val = {}
	various = {}
	
	# parse
	phase2Parse(input_stream, run,tags,runs,tag_runs,tags2,tag_val,various,1)
	# reexpoort

	if tsv_file == "":
		g = sys.stdout
		phase2Export(g,runs, tags,tag_runs,tags2,tag_val,nReadsGiven,nBasesGiven) 
	else:
		g = open(tsv_file+".aliqc.tsv", "w")
		phase2Export(g,runs, tags,tag_runs,tags2,tag_val,nReadsGiven,nBasesGiven) 
		g.close()

	sys.exit(0)
	
	
###############################################################################################	
###############################################################################################
# Phase3 parse the output and reexport in a nice format
###############################################################################################	

def viewHtmlMetaInfo (f) :
	f.write ( "  <META NAME=\"title\"\n")
	f.write (    " CONTENT=\"\n")
	f.write (  "Alignment statistics in next-generation sequencing\">\n\n")

	f.write ( "<META NAME=\"keywords\"\n")
	f.write (  " CONTENT=\"\n")
	f.write ( "Quality control BAM SAM")
	f.write (  "\">\n" )
	
	# this text is duplicated in index.html 
	f.write ( "<META NAME=\"Description\" \n")
	f.write ( " CONTENT= \"\n")
	f.write ( "aliqc.py is an open source script to generate alignemnt statistics from BAM/SAM files")
	f.write ( "\">\n")
	f.write ( "<meta http-equiv=\"content-type\" content=\"text/html;charset=iso-8859-1\">\n")

	f.write ( "<meta name=\"author\"\n")
	f.write ( " content=\"Joe Meehan, Danielle Thierry-Mieg and Jean Thierry-Mieg,\n")
	f.write ( "\">\n")

	# content privacy policy needed to be able to set a cookie 
	f.write ( "<meta http-equiv=\"P3P\" content='CP=\"IDC DSP COR CURa ADM OUR IND PHY ONL COM STA\"'>")

def viewOpenDocument (f, view_type, title):
	if view_type == "html" :
		f.write ("<html>\n")
		viewHtmlMetaInfo(f)
		f.write ("\n<body bgcolor=white>\n<h2>\n"+title+"\n</h2>\n")
	else:
		f.write ("# " + title + "\n")

def viewCloseDocument (f, view_type):
	if view_type == "html" :
		f.write ("\n</body>\n</html>\n") 
	else:
		f.write ("\n")

def viewTableOpen (f, view_type):
	if view_type == "html" :
		f.write ("<p>\n<table border=1>\n")
	else:
		f.write ("\n")

def viewTableClose (f, view_type,footnote):
	if view_type == "html" :
		f.write ("</table>\n"+footnote+"\n")
	else:
		f.write (footnote+"\n\n\n")

def viewTableNewLine (f, view_type, n):
	if view_type == "html" and n == 0:
		f.write ("<tr bgcolor=#afafff>\n") ;
	elif view_type == "html" and n % 2 == 0:
		f.write ("<tr bgcolor=#efefff>\n") ;
	elif view_type == "html" and n % 2 == 1:
		f.write ("<tr>\n") ;
	return n+1

def viewTableBreak (f, view_type):
	if view_type == "html" :
		f.write ("\n</tr>\n")
	else:
		f.write ("\n")
  
def viewTableCell (f, view_type, txt) :
	if view_type == "html" :
		f.write ("\t<td>"+txt+"</td>\n")
	else :
		f.write (tab+txt)

# if f == 0, we can use this call to get the numerical value of the tag
def viewTag (tag,f,run,tag_runs,tag_val,n) :
        z = "0"
        tag_run = tag + tab + run
        if tag_run in tag_runs and tag_runs[tag_run] > 0 :
                tag2 = tag_run + "#" + str(n)
                z = str(tag_val[tag2])
        if f :
                viewTableCell (f,view_type,z)
        return int(z)

def vBreak (f, view_type, txt) :
	if view_type == "html" :
		f.write ("<br>\n")
	else :
		f.write ("\n")

def phase3ExportStrand(f,pp,strand,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven ) :
 
        subkeys = ["AG","TC","GA","CT","AT","TA","GC","CG","AC","TG","GT","CA"]
        indelkeys = ["A", "T", "G", "C"]
	
        if (pp == 3 or pp == 4) and strand > 0:
                return

        if strand == 0 :
                viewTableOpen (f, view_type)

                if pp >= 3 and strand == 0 :
                        f.write("* In paired end sequencing, the position of the mismatches in the first and second read are given in succession. ")
                        f.write("This table may evidence cycles where the base call quality is poor and the mismatches are concentrated. ")
                        f.write("Althougth this is not apparent in the fastQ quality scores, this is frequently the case in Illumina sequencing, ")
                        f.write("and one could consider replacing the a-priori fastQ quality coefficient by an a-posteriori probability derived from the present table, as done by Magic, with good impact on SNP calling.\n")
                elif pp < 3 and strand == 0 :
                        f.write("# All reads, single or paired end.")
                elif strand == 1 :
                        f.write("# First reads in pairs.")
                elif strand == 2 :
                        f.write("# Second reads in pairs.")
		
                        if pp == 2:
                                f.write("### The mismatches are given in the orientation of the read.")


        n = viewTableNewLine (f, view_type, 0)

        if strand == 0 :
                viewTableCell (f, view_type, "# Read")
                viewTableCell (f, view_type, "Run")
                viewTableCell (f, view_type, "Method")
                if pp == 1 :
                        viewTableCell (f, view_type, "% aligned reads")
                        viewTableCell (f, view_type, "% aligned bases")
                        viewTableCell (f, view_type, "% matching bases")
                        viewTableCell (f, view_type, "% mapped pairs")
                        viewTableCell (f, view_type, "% compatible pairs*")
                        viewTableCell (f, view_type, "% fragments mapped on plus strand**")
                        viewTableCell (f, view_type, "% aligned reads with unique alignments")
                        viewTableCell (f, view_type, "% aligned reads with multiple alignments")
			
                        viewTableCell (f, view_type, "Reads in run or in file")
                        viewTableCell (f, view_type, "Reads unaligned")
                        viewTableCell (f, view_type, "Reads aligned")
			
                elif pp == 2 :
                        viewTableCell (f, view_type, "Reads in run or in file")
                        viewTableCell (f, view_type, "% Perfect reads exactly matching over their entire length")
                        viewTableCell (f, view_type, "% reads aligned with no mismatch")
                        viewTableCell (f, view_type, "% reads with 0 or 1 mismatch")

                if pp == 1 :
                        viewTableCell (f, view_type, "Reads with several alignments***")
                        imax= various["maxMultiAli"]
                        viewTableCell (f, view_type, "Reads with unique alignments")
                        for i in range(1,imax)  :
                                viewTableCell (f, view_type, "Reads with " +  str(i+1) + " alignments")
				
                        viewTableCell (f, view_type, "Fragments aligned on plus strand**")
                        viewTableCell (f, view_type, "Fragments aligned on minus strand**")
                        viewTableCell (f, view_type, "Pairs with both ends aligned")
                        viewTableCell (f, view_type, "Compatible pairs")
			
                        viewTableCell (f, view_type, "Raw bases in unaligned reads present in the file")
                        viewTableCell (f, view_type, "Raw bases in run or file")
		
                viewTableCell (f, view_type, "Aligned bases")

                if pp == 2 :
                        viewTableCell (f, view_type, "mismatches per kb aligned")
                        viewTableCell (f, view_type, "substitutions per kb aligned")
                        viewTableCell (f, view_type, "transitions per kb aligned")
                        viewTableCell (f, view_type, "transversions per kb aligned")
                        viewTableCell (f, view_type, "insertions per kb aligned")
                        viewTableCell (f, view_type, "deletions per kb aligned")
			
                        viewTableCell (f, view_type, "Perfect reads exactly matching over their entire length")
                        viewTableCell (f, view_type, "Reads with 0 mismatch")
                        imax= various["maxError"]
                        if imax > 0 :
                                viewTableCell (f, view_type, "1 mismatch")
                        for i in range(1,imax)  :
                                viewTableCell (f, view_type, str(i+1) + " mismatches")
                        viewTableCell (f, view_type, "Mismatch events")

                        viewTableCell (f, view_type, "Transitions counted in protected zones")
                        viewTableCell (f, view_type, "Transversions counted in protected zones")
                        viewTableCell (f, view_type, "Substitution events")
                        viewTableCell (f, view_type, "Single substitutions")
                        viewTableCell (f, view_type, "Double substitutions")
                        viewTableCell (f, view_type, "Triple substitutions")
                        viewTableCell (f, view_type, "Quadruple substitutions")
                        viewTableCell (f, view_type, "Quintuple substitutions")
                        viewTableCell (f, view_type, "Sextuple substitutions")
                        viewTableCell (f, view_type, "Septuple substitutions")
                        viewTableCell (f, view_type, "Octuple substitutions")
                        viewTableCell (f, view_type, "Longer 3 substitutions")
                        viewTableCell (f, view_type, "Inserted bases")
                        viewTableCell (f, view_type, "Insertion events")
                        viewTableCell (f, view_type, "Single insertion events")
                        viewTableCell (f, view_type, "Double insertion events")
                        viewTableCell (f, view_type, "Triple insertion events")
                        viewTableCell (f, view_type, "Quadruple insertion events")
                        viewTableCell (f, view_type, "Quintuple insertion events")
                        viewTableCell (f, view_type, "Sexuple insertion events")
                        viewTableCell (f, view_type, "Septuple insertion events")
                        viewTableCell (f, view_type, "Octuple insertion events")
                        viewTableCell (f, view_type, "Insertion events 9 to 40")
                        viewTableCell (f, view_type, "Insertion events 41 to 1000")
                        viewTableCell (f, view_type, "Insertion events over 1000")

                        viewTableCell (f, view_type, "Deleted base events")
                        viewTableCell (f, view_type, "Deletion event events")
                        viewTableCell (f, view_type, "Single deletion events")
                        viewTableCell (f, view_type, "Double deletion events")
                        viewTableCell (f, view_type, "Triple deletion events")
                        viewTableCell (f, view_type, "Quadruple deletion events")
                        viewTableCell (f, view_type, "Quintuple deletion events")
                        viewTableCell (f, view_type, "Sexuple deletion events")
                        viewTableCell (f, view_type, "Septuple deletion events")
                        viewTableCell (f, view_type, "Octuple deletion events")
                        viewTableCell (f, view_type, "Deletion events 9 to 40")
                        viewTableCell (f, view_type, "Deletion events 41 to 1000")
                        viewTableCell (f, view_type, "Deletion events over 1000")
                        viewTableCell (f, view_type, "Introns (usually gt-ag)")

                        viewTableCell (f, view_type, "Insertion events in polymers")
                        viewTableCell (f, view_type, "Single insertion events in polymers")
                        viewTableCell (f, view_type, "Double insertion events in polymers")
                        viewTableCell (f, view_type, "Triple insertion events in polymers")
                        viewTableCell (f, view_type, "Quadruple insertion events in polymers")
                        viewTableCell (f, view_type, "Quintuple insertion events in polymers")
                        viewTableCell (f, view_type, "Sexuple insertion events in polymers")
                        viewTableCell (f, view_type, "Septuple insertion events in polymers")
                        viewTableCell (f, view_type, "Octuple insertion events in polymers")
 
                        viewTableCell (f, view_type, "Deletion events in polymers")
                        viewTableCell (f, view_type, "Single deletion events in polymers")
                        viewTableCell (f, view_type, "Double deletion events in polymers")
                        viewTableCell (f, view_type, "Triple deletion events in polymers")
                        viewTableCell (f, view_type, "Quadruple deletion events in polymers")
                        viewTableCell (f, view_type, "Quintuple deletion events in polymers")
                        viewTableCell (f, view_type, "Sexuple deletion events in polymers")
                        viewTableCell (f, view_type, "Septuple deletion events in polymers")
                        viewTableCell (f, view_type, "Octuple deletion events in polymers")

                        viewTableCell (f, view_type, "Soft clipped bases")
                        viewTableCell (f, view_type, "Soft clipped events")
                        viewTableCell (f, view_type, "Hard clipped bases")
                        viewTableCell (f, view_type, "Hard clipped events")

                        viewTableCell (f, view_type, "Short 5prime overhang (1 to 16 bases)")
                        viewTableCell (f, view_type, "Mid 5prime overhang (17 to 100 bases)")
                        viewTableCell (f, view_type, "Long 5prime overhang (over 100 bases)")
                        viewTableCell (f, view_type, "Short 3prime overhang (1 to 16 bases)")
                        viewTableCell (f, view_type, "Mid 3prime overhang (17 to 100 bases)")
                        viewTableCell (f, view_type, "Long 3prime overhang (over 100 bases)")


                        for k in subkeys:
                                viewTableCell (f,view_type, k[0]+">"+k[1])
                        for k in indelkeys:
                                viewTableCell (f,view_type, "Ins "+k)
                        for k in indelkeys:
                                viewTableCell (f,view_type, "Del "+k)

                if pp == 1 :
                        imin= various["minAliLn"]
                        imax= various["maxAliLn"]+1
                        if imin < imax :
                                viewTableCell (f, view_type, "Aligned length")
                                for i in range(imin,imax) :
                                        viewTableCell (f, view_type, str(imax + imin - i -1)+" nt")

                if pp >= 3 :
                        imax= various["maxErrPos"] ;
                        imax1= various["maxErrPos1"]
                        imax2= various["maxErrPos2"]
                        if imax > 0 :
                                aliB1 = viewTableCell (f, view_type, "Aligned bases in read 1")
                                aliB2 = viewTableCell (f, view_type, "Aligned bases in read 2")
                        if pp == 3 :
                                viewTableCell (f, view_type, "Mismatch per sequencing cycle in read 1")
                        if pp == 4 :
                                viewTableCell (f, view_type, "Rate of Mismatch per Mb per sequencing cycle in read 1")
                        for i in range(imax+5)  :
                                viewTableCell (f, view_type, "Cycle " + str(i+1))
                        if imax2 > 0 :
                                if pp == 3 :
                                        viewTableCell (f, view_type, "Mismatch per sequencing cycle in read 2")
                                if pp == 4 :
                                        viewTableCell (f, view_type, "Rate of Mismatch per Mb per sequencing cycle in read 2")
                                for i in range(imax2+5)  :
                                        viewTableCell (f, view_type, "Cycle " + str(i+1))

        for run in sorted(runs):
                if strand > 0 : 
                        tag2 = "Reads_in_file" + tab + run + "#3"
                        nrf = 0
                        if tag2 in tag_val :
                                nrf = tag_val[tag2]
                                if nrf == 0 :
                                        continue

                viewTableBreak (f, view_type)
                n = viewTableNewLine (f, view_type, n)

                if strand == 0 :
                        viewTableCell (f, view_type, "Any")
                elif strand == 1 :
                        viewTableCell (f, view_type, "R1")
                elif strand == 2 :
                        viewTableCell (f, view_type, "R2")

                runBuf = run.split(".")
                ndot= len(runBuf)
                if (ndot == 1) :
                        viewTableCell (f, view_type, runBuf[0])
                        viewTableCell (f, view_type, "0")
                if (ndot == 2) :
                        viewTableCell (f, view_type, runBuf[0] + "." + runBuf[1])
                        viewTableCell (f, view_type, "0")
                if (ndot == 3) :
                        viewTableCell (f, view_type, runBuf[0] + "." + runBuf[1])
                        viewTableCell (f, view_type, runBuf[2])
                if (ndot >= 4) :
                        viewTableCell (f, view_type, runBuf[0] + "." + runBuf[1])
                        viewTableCell (f, view_type, runBuf[2] + "." + runBuf[3])

                if pp == 1 :
                        r1 = viewTag ("Reads_in_file",0,run,tag_runs,tag_val,1+strand)
                        r1a = viewTag ("Reads_in_run",0,run,tag_runs,tag_val,1+strand)
                        if nReadsGiven :
                                r1a = nReadsGiven
                                if strand > 0 :
                                        r1a /= 2 
                        if r1a > 0 :
                                r1 = r1a

                                r2 = viewTag ("Reads_Mapped",0,run,tag_runs,tag_val,1+strand)
                        if r1 == 0 :
                                r1 = 1
                        r = 100.0 * r2/r1
                        viewTableCell (f, view_type, "%.2f" % (r))
			
                        r3 = viewTag ("Cumulated_read_length",0,run,tag_runs,tag_val,1+strand)
                        if r2 > 0 and r1 > r2:
                                r3 = r3 * r1/r2
                        r3a = viewTag ("Bases_in_run",0,run,tag_runs,tag_val,1+strand)
                        if nBasesGiven > 0:
                                r3a = nBasesGiven
                                if strand > 0 :
                                        r3a /= 2 
                        if r3a > 0 :
                                r3 = r3a
                        if r3 == 0 :
                                r3 = 1

                        r4 = viewTag ("Aligned_bases",0,run,tag_runs,tag_val,1+strand)
                        r4b = viewTag ("Matched_bases",0,run,tag_runs,tag_val,1+strand)
                        r4c = viewTag ("Michmatched_bases",0,run,tag_runs,tag_val,1+strand)
                        if r4c > 0:
                                r4 = r4b + r4c
                        r = 100.0 * r4/r3
                        viewTableCell (f, view_type, "%.2f" % (r))
                        r = 100.0 * r4b/r3
                        viewTableCell (f, view_type, "%.2f" % (r))


                        sp0 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 1 + 2 * 0)
                        sp1 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 1 + 2 * 1)
                        sp2 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 1 + 2 * 2)
                        sm0 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 2 + 2 * 0)
                        sm1 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 2 + 2 * 1)
                        sm2 = viewTag ("Reads_Mapped_per_strand",0,run,tag_runs,tag_val, 2 + 2 * 2)

                        if sp1 + sp2 == 0 :
                                sp1 = sp0
                        if sm1 + sm2 == 0 :
                                sm1 = sm0

                        if strand == 0 :
                                sx = sp1 + sm2
                                sy = sm1 + sp2
                        elif strand == 1 :
                                sx = sp1
                                sy = sm1
                        elif strand == 2 :
                                sx = sm2
                                sy = sp2

                        if strand == 0 :
                                r5 = viewTag ("Both_Reads_Mapped",0,run,tag_runs,tag_val,1)/2
                                r = 200.0 * r5/(r1)
                                viewTableCell (f, view_type, "%.2f" % (r))
                                r6 = viewTag ("Both_Reads_Mapped",0,run,tag_runs,tag_val,2)/2
                                r = 200.0 * r6/(r1)
                                viewTableCell (f, view_type, "%.2f" % (r))
                        else:
                                viewTableCell (f, view_type, "-")
                                viewTableCell (f, view_type, "-")

                        if sx == 0 :
                                r = 0
                        else :
                                r = 100.0 * sx/(sx + sy)
                        viewTableCell (f, view_type, "%.4f" % (r))

                        r7 = viewTag (("MultiAli:%03d" % (1)),0,run,tag_runs,tag_val,1+strand) 
                        if r2 == 0 :
                                ua = 0
                                una = 0
                        else :
                                ua = 100.0 * r7/r2
                                una = 100 - ua

                        viewTableCell (f, view_type, "%.2f" % (ua))
                        viewTableCell (f, view_type, "%.2f" % (una))

                        viewTableCell (f, view_type, "%d" % (r1))
                        r = r1 - r2  # r1: reads in run, r2: reads aligneed, r1-r2: reads unaligned
                        viewTableCell (f, view_type, "%d" % (r))	
                        viewTableCell (f, view_type, "%d" % (r2))

                        r = r2 - r7
                        viewTableCell (f, view_type, "%d" % (r))   # reads with several ali	
                        # viewTag ("Alignments",f,run,tag_runs,tag_val,1+strand)
                        viewTableCell (f, view_type, "%d" % (r7))	
                        imax= various["maxMultiAli"]
                        for i in range(1,imax) :
                                viewTag (("MultiAli:%03d" % (i+1)),f,run,tag_runs,tag_val,1+strand) 

                        viewTableCell (f, view_type, "%d" % (sx))
                        viewTableCell (f, view_type, "%d" % (sy))
			
                        if strand == 0 :
                                r = int(viewTag ("Both_Reads_Mapped",0,run,tag_runs,tag_val,1)/2)
                                viewTableCell (f, view_type, str(r))
                                r = int (viewTag ("Both_Reads_Mapped",0,run,tag_runs,tag_val,2)/2)
                                viewTableCell (f, view_type, str(r))
                        else:
                                viewTableCell (f, view_type, "-")
                                viewTableCell (f, view_type, "-")
                        r = viewTag ("Unaligned_bases",0,run,tag_runs,tag_val,1+strand)
                        if r == 0 :
                                viewTableCell (f, view_type, "-")
                        else:
                                viewTag ("Unaligned_bases",f,run,tag_runs,tag_val,1+strand)
                        viewTableCell (f, view_type, "%d" % (r3))   # nBasesGiven raw bases in run or file

                if pp == 2 :
                        zr = viewTag ("Reads_in_file",0,run,tag_runs,tag_val,1+strand)
                        zr2 = viewTag ("Reads_in_run",0,run,tag_runs,tag_val,1+strand)
                        if nReadsGiven :
                                zr2 = nReadsGiven
                                if strand > 0 :
                                        zr2 /= 2 
                        if zr2 > 0 :
                                zr = zr2
                        viewTableCell (f, view_type, "%d" % (zr))
                        if zr == 0 :
                                zr = 1
                        z = viewTag ("Perfect_reads",0,run,tag_runs,tag_val,1+strand)
                        r = 100.0 * z/zr
                        viewTableCell (f, view_type, "%.2f" % (r))
                        z0 = viewTag (("MismatchPerRead:%03d" % (0)),0,run,tag_runs,tag_val,1+strand)
                        r = 100.0 * z0/zr
                        viewTableCell (f, view_type, "%.2f" % (r))
                        z1 = viewTag (("MismatchPerRead:%03d" % (1)),0,run,tag_runs,tag_val,1+strand)
                        r = 100.0 * (z0 + z1)/zr
                        viewTableCell (f, view_type, "%.2f" % (r))

                viewTag ("Aligned_bases",f,run,tag_runs,tag_val,1+strand)

                if pp == 1 :
                        imin= various["minAliLn"]
                        imax= various["maxAliLn"]+1
                        if imin < imax :
                                viewTableCell (f, view_type, run)
                                for i in range(imin,imax) :
                                        z0 = viewTag (("AliLn:%04d" % (imax+imin-i-1)),0,run,tag_runs,tag_val,1) 
                                        z1 = viewTag (("AliLn:%04d" % (imax+imin-i-1)),0,run,tag_runs,tag_val,2) 
                                        z2 = viewTag (("AliLn:%04d" % (imax+imin-i-1)),0,run,tag_runs,tag_val,3) 
                                        if strand == 1 and z0 > z1 + z2 :
                                                viewTag (("AliLn:%04d" % (imax+imin-i-1)),f,run,tag_runs,tag_val,1) 
                                        else:
                                                viewTag (("AliLn:%04d" % (imax+imin-i-1)),f,run,tag_runs,tag_val,strand+1) 

                if pp == 2 :
                        strand2 = strand            
                        z = viewTag ("Alignments",0,run,tag_runs,tag_val,1+0)
                        z1 = viewTag ("Alignments",0,run,tag_runs,tag_val,1+1)
                        z2 = viewTag ("Alignments",0,run,tag_runs,tag_val,1+2)
                        if z1 + z2 == 0 and strand == 1 :
                                strand2 = 0
                        zab = viewTag ("Aligned_bases",0,run,tag_runs,tag_val,1+strand2)
                        if zab == 0 :
                                zab = 1
                        z = viewTag ("Mismatch_events",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))
                        z = viewTag ("Substitutions",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))
                        z = viewTag ("Transitions",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))
                        z = viewTag ("Transversions",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))
                        z = viewTag ("Inserted_bases",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))
                        z = viewTag ("Deleted_bases",0,run,tag_runs,tag_val,1+strand2)
                        r = 1000.0 * z/zab
                        viewTableCell (f, view_type, "%.3f" % (r))

                        viewTag ("Perfect_reads",f,run,tag_runs,tag_val,1+strand2)
                        imax= various["maxError"]
                        for i in range(imax+1)  :
                                viewTag (("MismatchPerRead:%03d" % (i)),f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Mismatch_events",f,run,tag_runs,tag_val,1+strand2)

                        viewTag ("Transitions",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Transversions",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Substitution_events",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Single_substitutions",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Double_substitutions",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Triple_substitutions",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Longer_substitutions",f,run,tag_runs,tag_val,1+strand2)



                        viewTag ("Inserted_bases",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Insertion_events",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Single_insertion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Double_insertion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Triple_insertion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Insertion_up_to_40",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Insertion_up_to_1000",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Insertion_long",f,run,tag_runs,tag_val,1+strand2)


                        viewTag ("Deleted_bases",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Deletion_events",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Single_deletion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Double_deletion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Triple_deletion",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Deletion_up_to_40",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Deletion_up_to_1000",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Deletion_long",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Intron",f,run,tag_runs,tag_val,1+strand2)

                        viewTag ("Insertion_events_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Single_insertion_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Double_insertion_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Triple_insertion_in_polymers",f,run,tag_runs,tag_val,1+strand2)

                        viewTag ("Deletion_events_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Single_deletion_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Double_deletion_in_polymers",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Triple_deletion_in_polymers",f,run,tag_runs,tag_val,1+strand2)

                        viewTag ("Hard_clipped_bases",f,run,tag_runs,tag_val,1)
                        viewTag ("Hard_clipped_events",f,run,tag_runs,tag_val,1)
                        viewTag ("Soft_clipped_bases",f,run,tag_runs,tag_val,1)
                        viewTag ("Soft_clipped_events",f,run,tag_runs,tag_val,1)

                        viewTag ("Unaligned_prefix_1_16",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Unaligned_prefix_17_100",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Unaligned_prefix_long",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Unaligned_suffix_1_16",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Unaligned_suffix_17_100",f,run,tag_runs,tag_val,1+strand2)
                        viewTag ("Unaligned_suffix_long",f,run,tag_runs,tag_val,1+strand2)

                        for k in subkeys:
                                viewTag (("Sub:%s" % (k)),f,run,tag_runs,tag_val,1+strand2)
                        for k in indelkeys:
                                viewTag (("Ins:%s" % (k)),f,run,tag_runs,tag_val,1+strand2)
                        for k in indelkeys:
                                viewTag (("Del:%s" % (k)),f,run,tag_runs,tag_val,1+strand2)




                if pp >= 3 :
                        imax = various["maxErrPos"]
                        imax1 = various["maxErrPos1"]
                        imax2 = various["maxErrPos2"]
                        aliB = 1
                        aliB1 = 1
                        aliB2 = 1
                        hack = 0
                        if imax > 0 :
                                aliB = viewTag ("Aligned_bases",0,run,tag_runs,tag_val,1+0)
                                aliB1 = viewTag ("Aligned_bases",0,run,tag_runs,tag_val,1+1)
                                aliB2 = viewTag ("Aligned_bases",0,run,tag_runs,tag_val,1+2)
                        if (aliB1 + aliB2 == 0) :
                                hack=1
                                aliB = viewTag ("Aligned_bases",f,run,tag_runs,tag_val,1+0)
                        else:
                                aliB1 = viewTag ("Aligned_bases",f,run,tag_runs,tag_val,1+1)
                        aliB2 = viewTag ("Aligned_bases",f,run,tag_runs,tag_val,1+2)

                        if (aliB == 0) :
                                aliB = 1
                        if (aliB1 + aliB2 == 0) :
                                aliB1 = aliB
                        if (aliB2 == 0) :
                                aliB2 = 1
                        viewTableCell (f, view_type, "")
                        for i in range(imax + 5) :
                                z = viewTag (("Err_pos:%04d" % (i)),0,run,tag_runs,tag_val,1+0)
                                z1 = viewTag (("Err_pos:%04d" % (i)),0,run,tag_runs,tag_val,1+1)
                                z2 = viewTag (("Err_pos:%04d" % (i)),0,run,tag_runs,tag_val,1+2)
                                if z1 + z2 > 0 :
                                        z = z1
                                if (pp == 3):
                                        viewTableCell (f,view_type, ("%d" % (z)))
                                if (pp == 4):
                                        viewTableCell (f,view_type, ("%.2f" % (1000000.0*z/aliB1)))
                        if imax2 > 0 :
                                viewTableCell (f, view_type, "")
                                for i in range(imax2 + 5) :
                                        z = viewTag (("Err_pos:%04d" % (i)),0,run,tag_runs,tag_val,1+2)
                                        if (pp == 3):
                                                viewTableCell (f,view_type, ("%d" % (z)))
                                        if (pp == 4):
                                                viewTableCell (f,view_type, ("%.2f" % (1000000.0*z/aliB2)))

        viewTableBreak (f, view_type)

        if strand == 2 :
                footnote = "### The definition of a compatible pair depends on the aligner, it has no absolute meaning"	
                if view_type == "html" :
                        footnote = footnote + "<br>"
                footnote = footnote + "\n### In paired end sequencing, the strand specificity of the fragments is evaluated after complementing the R2 read, since the R1 and R2 reads face each other"
                if view_type == "html" :
                        footnote = footnote + "<br>"
                footnote = footnote + "\n### A read may have several alignments, but each read contributes only once to the other columns of this table, on its first occurence in the alignment file"
                viewTableClose (f, view_type,footnote)
	

def phase3date (f, title) :
	txt = time.strftime ("%Y-%m-%d %H:%M:%S")
	if view_type == "html" :
		f.write ("#<i> "+txt+"</i>+"     "+title<br>\n")
	else :
		f.write ("# "+txt+"\t"+title+"\n")
				

def phase3Export(f,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven,multi_file) :

        if multi_file == 0 :
                viewOpenDocument (f, view_type, "Alignement statistics, there are 4 tables in this page")
                phase3date (f, "Alignement statistics,  minimal alignment " + str(minAli) + " bp") 

        if multi_file > 0 :
                f = open(aligned_reads_multiplicity_length_aligned_file,"w")
                phase3date (f, "Aligned reads multiplicity length aligned,  minimal alignment " + str(minAli) + " bp") 
        for strand in range(3) :
                phase3ExportStrand(f,1,strand,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven)

        if multi_file > 0 :
                f.close()
                f = open(mismatch_histo_and_types_file,"w")
                phase3date (f, "Mismatch histo and types,  minimal alignment " + str(minAli) + " bp")
        for strand in range(3) :
                phase3ExportStrand(f,2,strand,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven)
	
        if 0:
                if multi_file > 0 :
                        f.close()
                        f = open(mismatch_counts_per_cycle_file,"w")
                        phase3date (f, "Mismatch counts per cycle,  minimal alignment " + str(minAli) + " bp") 
                for strand in range(1) :
                        phase3ExportStrand(f,3,strand,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven)

        if multi_file > 0 :
                f.close()
                f = open(mismatch_per_cycle_per_kb_aligned_file,"w")
                phase3date (f, "Mismatch per cycle per kb aligned,  minimal alignment " + str(minAli) + " bp") 
        for strand in range(1) :
                phase3ExportStrand(f,4,strand,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven)

        viewCloseDocument (f, view_type)
	
###############################################################################################

def phase3(input_stream, view_type, out_prefix, nReadsGiven,nBasesGiven,multi_file):
        tags = {}
        tags2 = {}
        runs = {}
        tag_val = {}
        tag_runs = {}
        tag_val = {}
        various = {}

        # parse
        phase2Parse(input_stream, "any",tags,runs,tag_runs,tags2,tag_val,various,0)
        # parse the input, recognize the runId in column 1. the tags in col 2, add the counts
        # reexpoort
        g = 0
        if out_prefix == "":
                g = sys.stdout
                phase3Export(g,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven,multi_file)
        else:
                if multi_file == 0:
                        if view_type == "html":
                                g = open(out_prefix+".html", "w")
                        if view_type == "table":
                                g = open(out_prefix+".tsv", "w")
        phase3Export(g,view_type,runs, tags,tag_runs,tags2,tag_val,various,nReadsGiven,nBasesGiven,multi_file)
        if multi_file == 0:
                g.close()

        sys.exit(0)
	
	
###############################################################################################	
###############################################################################################	
# exportation functions
###############################################################################################	
# summary

def exportSummary (f,run_id, view_type, nReadsGiven,nBasesGiven) :

        global fasta_file
        global record_count
        global reads_in_file
        global reads_1_in_file
        global reads_2_in_file
        global alignment_count
        global perfect
        global perfect1
        global perfect2
        global perfect_clipped
        global perfect_clipped1
        global perfect_clipped2
        global read_bp
        global read_1_bp
        global read_2_bp
        global readU_bp
        global readU_1_bp
        global readU_2_bp
        global ali_bp
        global ali_1_bp
        global ali_2_bp
        global match_base
        global sub_base
        global del_base
        global ins_base
        global match_base_1
        global sub_base_1
        global del_base_1
        global ins_base_1
        global match_base_2
        global sub_base_2
        global del_base_2
        global ins_base_2
        global read1_alignment_count
        global read2_alignment_count
        global unaligned_count
        global read1_unaligned_count
        global read2_unaligned_count
	
        # Calculate transitions from base matrix totals
        transition_list = ["AG","GA","CT","TC"]
        transitions = 0
        transitions1 = 0
        transitions2 = 0
        for t in transition_list:
                if not t in base_matrix_total : # jfm 9/28/2012 total1 & total2 may not have all the keys total has
                        base_matrix_total[t] = 0
                if not t in base_matrix_total1 :
                        base_matrix_total1[t] = 0
                if not t in base_matrix_total2 :
                        base_matrix_total2[t] = 0
                transitions += base_matrix_total[t]
                transitions1 += base_matrix_total1[t]
                transitions2 += base_matrix_total2[t]
		
                # Calculate transversions from base matrix totals
        transversion_list = ["AC","CA","GT","TG","AT","TA","CG","GC"]
        transversions = 0
        transversions1 = 0
        transversions2 = 0
        for t in transversion_list:
                if not t in base_matrix_total : # jfm 9/28/2012 total1 & total2 may not have all the keys total has
                        base_matrix_total[t] = 0
                if not t in base_matrix_total1 :
                        base_matrix_total1[t] = 0
                if not t in base_matrix_total2 :
                        base_matrix_total2[t] = 0
                transversions += base_matrix_total[t]
                transversions1 += base_matrix_total1[t]
                transversions2 += base_matrix_total2[t]

        if view_type  == "tsv":
                f.write("Reference_File:"+fasta_file+tab+run_id+"\t0"+newline)
                f.write("Reads_in_file"+tab+run_id+"\t3"+tab+str(reads_in_file)+tab+str(reads_1_in_file)+tab+str(reads_2_in_file)+newline)

                nReadsGiven1 = 0
                nReadsGiven2 = 0
                if nReadsGiven > 0 :
                        nReadsGiven1 = nReadsGiven
                        if reads_1_in_file > 0 and reads_2_in_file > 0 :
                                nReadsGiven1 = nReadsGiven/2
                                nReadsGiven2 = nReadsGiven/2
                        elif reads_1_in_file > 0 and reads_2_in_file == 0 :
                                nReadsGiven1 = nReadsGiven
                                nReadsGiven2 = 0
                        elif reads_1_in_file == 0 and reads_2_in_file > 0 :
                                nReadsGiven2 = nReadsGiven
                                nReadsGiven1 = 0

                        f.write("Reads_in_run"+tab+run_id+"\t3"+tab+str(nReadsGiven)+tab+str(nReadsGiven1)+tab+str(nReadsGiven2)+newline)

                nBasesGiven1 = 0
                nBasesGiven2 = 0
                if nBasesGiven > 0 :
                        nBasesGiven1 = nBasesGiven
                        if reads_1_in_file > 0 and reads_2_in_file > 0 :
                                nBasesGiven1 = nBasesGiven/2
                                nBasesGiven2 = nBasesGiven/2
                        elif reads_1_in_file > 0 and reads_2_in_file == 0 :
                                nBasesGiven1 = nBasesGiven
                                nBasesGiven2 = 0
                        elif reads_1_in_file == 0 and reads_2_in_file > 0 :
                                nBasesGiven2 = nBasesGiven
                                nBasesGiven1 = 0

                        f.write("Bases_in_run"+tab+run_id+"\t3"+tab+str(nBasesGiven)+tab+str(nBasesGiven1)+tab+str(nBasesGiven2)+newline)

                f.write("Alignments"+tab+run_id+"\t3\t"+str(alignment_count)+tab+str(read1_alignment_count)+tab+str(read2_alignment_count)+newline)
                f.write("Unaligned_reads"+tab+run_id+"\t3\t"+str(unaligned_count)+tab+str(read1_unaligned_count)+tab+str(read2_unaligned_count)+newline)
                f.write("Cumulated_read_length"+tab+run_id+"\t3\t"+str(read_bp)+tab+str(read_1_bp)+tab+str(read_2_bp)+newline)
                mismatch_base = sub_base + del_base + ins_base
                mismatch_base_1 = sub_base_1 + del_base_1 + ins_base_1
                mismatch_base_2 = sub_base_2 + del_base_2 + ins_base_2
                ali_bp = match_base + mismatch_base
                ali_1_bp = match_base_1 + mismatch_base_2
                ali_2_bp = match_base_2 + mismatch_base_2
                f.write("Aligned_bases"+tab+run_id+"\t3\t"+str(ali_bp)+tab+str(ali_1_bp)+tab+str(ali_2_bp)+newline)
                f.write("Matched_bases"+tab+run_id+"\t3\t"+str(match_base)+tab+str(match_base_1)+tab+str(match_base_2)+newline)
                f.write("Misatched_bases"+tab+run_id+"\t3\t"+str(mismatch_base)+tab+str(mismatch_base_1)+tab+str(mismatch_base_2)+newline)
                f.write("Substitutions"+tab+run_id+"\t3\t"+str(sub_base)+tab+str(sub_base_1)+tab+str(sub_base_2)+newline)
                f.write("Transitions"+tab+run_id+"\t3\t"+str(transitions)+tab+str(transitions1)+tab+str(transitions2)+newline)		
                f.write("Transversions"+tab+run_id+"\t3\t"+str(transversions)+tab+str(transversions1)+tab+str(transversions2)+newline)
                f.write("Deleted_bases"+tab+run_id+"\t3\t"+str(del_base)+tab+str(del_base_1)+tab+str(del_base_2)+newline)
                f.write("Inserted_bases"+tab+run_id+"\t3\t"+str(ins_base)+tab+str(ins_base_1)+tab+str(ins_base_2)+newline)
                f.write("Unaligned_bases"+tab+run_id+"\t3\t"+str(readU_bp)+tab+str(readU_1_bp)+tab+str(readU_2_bp)+newline)
		
        global reads_mapped_count
        global read1_mapped_count
        global read2_mapped_count
        global pair_mapped_count
        global proper_pair_count

        if view_type  == "tsv":
                f.write("Reads_Mapped"+tab+run_id+"\t3\t"+str(0+reads_mapped_count)+tab+str(0+read1_mapped_count)+tab+str(0+read2_mapped_count)+newline)
                f.write("Both_Reads_Mapped"+tab+run_id+"\t2\t"+str(pair_mapped_count)+tab+str(proper_pair_count)+newline)
        global reads_plus_strand_count
        global reads_minus_strand_count
        global read1_plus_strand_count
        global read1_minus_strand_count
        global read2_plus_strand_count
        global read2_minus_strand_count

        if view_type  == "tsv":
                f.write("Reads_Mapped_per_strand"+tab+run_id+"\t6")
                f.write(tab+str(reads_plus_strand_count))
                f.write(tab+str(reads_minus_strand_count))
                f.write(tab+str(read1_plus_strand_count))
                f.write(tab+str(read1_minus_strand_count))
                f.write(tab+str(read2_plus_strand_count))
                f.write(tab+str(read2_minus_strand_count)+newline)
                f.write("Perfect_reads"+tab+run_id+"\t3"+tab+str(perfect)+tab+str(perfect1)+tab+str(perfect2)+newline)
                f.write("Perfect_clipped"+tab+run_id+"\t3"+tab+str(perfect_clipped)+tab+str(perfect_clipped1)+tab+str(perfect_clipped2)+newline)
		
###############################################################################################	
# aligned length histogram

def exportAlignedLngths(f,run_id,view_type,match_length_hist) :
	keys = sorted(list(match_length_hist.keys()))
	for k in keys:
		if not k in match_length_hist1: # jfm 9/28/2012 hist1 & hist2 may not have all the keys hist has
			match_length_hist1[k] = 0
		if not k in match_length_hist2:
	       		match_length_hist2[k] = 0
		if view_type  == "tsv":
			f.write (("AliLn:%04d" % (k))+tab+run_id+"\t3\t"+str(match_length_hist[k])+tab+str(match_length_hist1[k])+tab+str(match_length_hist2[k])+newline )
			
###############################################################################################	
# mismatch rate histogram

def exportAnyEvent(f,run_id,view_type,any_event_hist) :

        keys = sorted(list(any_event_hist.keys()))
        for k in keys:
                if not k in any_event_hist1: # jfm 9/28/2012 hist1 & hist2 may not have all the keys hist has
                        any_event_hist1[k] = 0
                if not k in any_event_hist2:
                        any_event_hist2[k] = 0
                if view_type  == "tsv":
                        f.write(("MismatchPerRead:%03d" % (k))+tab+run_id+"\t3\t"+str( any_event_hist[k])+tab+str( any_event_hist1[k])+tab+str( any_event_hist2[k])+newline )
	
###############################################################################################	
# error position histogram

def exportMismatchPositions(f,run_id,view_type,err_pos,err_pos1,err_pos2,err_N1,err_N2) :

        keys = sorted(set(list(err_pos.keys())+list(err_N1.keys())+list(err_N2.keys()))) # jfm 10/4/2012 added N1 and N2 keys
        for k in keys:
                n0 = 0
                n1 = 0
                n2 = 0
                nn0 = 0
                nn1 = 0
                nn2 = 0
                if k in err_pos :
                        n0 = err_pos[k]
                if k in err_pos1 :
                        n1 = err_pos1[k]
                if k in err_pos2 :
                        n2 = err_pos2[k]
                if k in err_N1 :
                        nn1 = err_N1[k]
                if k in err_N2 :
                        nn2 = err_N2[k]
                if n1 + n2 == 0:
                        n1 = n0
                if (n0 + n1 + n2 + nn1 + nn2 > 0) :
                        if view_type  == "tsv":
                                f.write(("Err_pos:%04d" %(k))+tab+run_id+"\t3\t"+str(n1+n2)+tab+str(n1)+tab+str(n2)+newline)
                                f.write(("Ambiguous_pos:%04d" %(k))+tab+run_id+"\t3\t"+str(nn1+nn2)+tab+str(nn1)+tab+str(nn2)+newline) 
	
###############################################################################################	
# insertions

def exportInsertion(f,run_id,  view_type,  insert_type_hist, insert_type_hist1, insert_type_hist2) :
        bases = ["A","T","G","C", "N"]
        keys = list( set (bases + list(insert_type_hist.keys())) )  # Make sure the single bases are in the list
        # keys.sort()
        for k in keys:
                if not k in insert_type_hist : # jfm 9/28/2012 some types may not occur
                        insert_type_hist[k] = 0
                if not k in insert_type_hist1 : 
                        insert_type_hist1[k] = 0
                if not k in insert_type_hist2 :
                        insert_type_hist2[k] = 0
                if view_type  == "tsv":
                        f.write("Ins:"+k+tab+run_id+"\t3\t"+str(insert_type_hist[k])+tab+str(insert_type_hist1[k])+tab+str(insert_type_hist2[k])+newline)
			
###############################################################################################	
# deletions

def exportDeletions(f,run_id,view_type,delete_type_hist, delete_type_hist1, delete_type_hist2 ):
	bases = ["A","T","G","C","N"]
	keys = list( set (bases + list(delete_type_hist.keys())) )  # Make sure the single bases are in the list
	# keys.sort()
	for k in keys:
		if not k in delete_type_hist : # jfm 9/28/2012 some types may not occur
			delete_type_hist[k] = 0
		if not k in delete_type_hist1 :
			delete_type_hist1[k] = 0
		if not k in delete_type_hist2 :
			delete_type_hist2[k] = 0
		if view_type  == "tsv":
			f.write("Del:"+k+tab+run_id+"\t3\t"+str(delete_type_hist[k])+tab+str(delete_type_hist1[k])+tab+str(delete_type_hist2[k])+newline)
			
###############################################################################################
# substitutions

def exportSubstitutions(f,run_id,view_type,base_matrix_total,base_matrix_total1,base_matrix_total2) :

        global hard_clipped_events
        global hard_clipped_bases
        global soft_clipped_events
        global soft_clipped_bases

        global intron_N
        global intron_N_1
        global intron_N_2
        
        global intron_n
        global intron_n_1
        global intron_n_2
        
        global prefix_S
        global prefix_S_1
        global prefix_S_2
        
        global suffix_S
        global suffix_S_1
        global suffix_S_2
        
        global prefix_M
        global prefix_M_1
        global prefix_M_2
        
        global suffix_M
        global suffix_M_1
        global suffix_M_2
        
        global prefix_L
        global prefix_L_1
        global prefix_L_2
        
        global suffix_L
        global suffix_L_1
        global suffix_L_2

        global any_substitutions
        global any_substitutions1
        global any_substitutions2
        
        global any_insertions
        global any_insertions1
        global any_insertions2

        global insertions_40
        global insertions_40_1
        global insertions_40_2

        global insertions_1000
        global insertions_1000_1
        global insertions_1000_2

        global insertions_long
        global insertions_long1
        global insertions_long2

        global deletions_40
        global deletions_40_1
        global deletions_40_2

        global deletions_1000
        global deletions_1000_1
        global deletions_1000_2

        global deletions_long
        global deletions_long1
        global deletions_long2

        global sliding_insertions
        global sliding_insertions1
        global sliding_insertions2

        global any_deletions
        global any_deletions1
        global any_deletions2

        global sliding_deletions
        global sliding_deletions1
        global sliding_deletions2

        global singlet_substitutions
        global singlet_substitutions1
        global singlet_substitutions2

        global singlet_insertions
        global singlet_insertions1
        global singlet_insertions2

        global singlet_deletions
        global singlet_deletions1
        global singlet_deletions2

        global doublet_substitutions
        global doublet_substitutions1
        global doublet_substitutions2

        global doublet_insertions
        global doublet_insertions1
        global doublet_insertions2

        global doublet_deletions
        global doublet_deletions1
        global doublet_deletions2

        global triplet_substitutions
        global triplet_substitutions1
        global triplet_substitutions2

        global longer_substitutions
        global longer_substitutions1
        global longer_substitutions2

        global triplet_insertions
        global triplet_insertions1
        global triplet_insertions2

        global triplet_deletions
        global triplet_deletions1
        global triplet_deletions2

        global quadruplet_substitutions
        global quadruplet_substitutions1
        global quadruplet_substitutions2
        
        global quadruplet_insertions
        global quadruplet_insertions1
        global quadruplet_insertions2
        
        global quadruplet_deletions
        global quadruplet_deletions1
        global quadruplet_deletions2
        
        global quintuplet_substitutions
        global quintuplet_substitutions1
        global quintuplet_substitutions2
        
        global quintuplet_insertions
        global quintuplet_insertions1
        global quintuplet_insertions2
        
        global quintuplet_deletions
        global quintuplet_deletions1
        global quintuplet_deletions2
        
        global sextuplet_substitutions
        global sextuplet_substitutions1
        global sextuplet_substitutions2
        
        global sextuplet_insertions
        global sextuplet_insertions1
        global sextuplet_insertions2
        
        global sextuplet_deletions
        global sextuplet_deletions1
        global sextuplet_deletions2
        
        global septuplet_substitutions
        global septuplet_substitutions1
        global septuplet_substitutions2
        
        global septuplet_insertions
        global septuplet_insertions1
        global septuplet_insertions2
        
        global septuplet_deletions
        global septuplet_deletions1
        global septuplet_deletions2
        
        global octuplet_substitutions
        global octuplet_substitutions1
        global octuplet_substitutions2
        
        global octuplet_insertions
        global octuplet_insertions1
        global octuplet_insertions2
        
        global octuplet_deletions
        global octuplet_deletions1
        global octuplet_deletions2

        global singlet_s_insertions
        global singlet_s_insertions1
        global singlet_s_insertions2

        global singlet_s_deletions
        global singlet_s_deletions1
        global singlet_s_deletions2

        global doublet_s_insertions
        global doublet_s_insertions1
        global doublet_s_insertions2

        global doublet_s_deletions
        global doublet_s_deletions1
        global doublet_s_deletions2

        global triplet_s_insertions
        global triplet_s_insertions1
        global triplet_s_insertions2

        global triplet_s_deletions
        global triplet_s_deletions1
        global triplet_s_deletions2

        global quadruplet_s_insertions
        global quadruplet_s_insertions1
        global quadruplet_s_insertions2
        
        global quadruplet_s_deletions
        global quadruplet_s_deletions1
        global quadruplet_s_deletions2
        
        global quintuplet_s_insertions
        global quintuplet_s_insertions1
        global quintuplet_s_insertions2
        
        global quintuplet_s_deletions
        global quintuplet_s_deletions1
        global quintuplet_s_deletions2
        
        global sextuplet_s_insertions
        global sextuplet_s_insertions1
        global sextuplet_s_insertions2
        
        global sextuplet_s_deletions
        global sextuplet_s_deletions1
        global sextuplet_s_deletions2
        
        global septuplet_s_insertions
        global septuplet_s_insertions1
        global septuplet_s_insertions2
        
        global septuplet_s_deletions
        global septuplet_s_deletions1
        global septuplet_s_deletions2
        
        global octuplet_s_insertions
        global octuplet_s_insertions1
        global octuplet_s_insertions2
        
        global octuplet_s_deletions
        global octuplet_s_deletions1
        global octuplet_s_deletions2
        
        global longer_s_insertions
        global longer_s_insertions1
        global longer_s_insertions2

        global longer_s_deletions
        global longer_s_deletions1
        global longer_s_deletions2

        global every_base_flag

	# Print (transitions, transversions, and base matrix totals)
        if view_type  == "tsv":
                f.write("Mismatch_events"+tab+run_id+"\t3\t"+str(any_substitutions+any_insertions+any_deletions)+tab+str(any_substitutions1+any_insertions1+any_deletions1)+tab+str(any_substitutions2+any_insertions2+any_deletions2)+newline)

                f.write("Substitution_events"+tab+run_id+"\t3\t"+str(any_substitutions)+tab+str(any_substitutions1)+tab+str(any_substitutions2)+newline)
                f.write("Single_substitutions"+tab+run_id+"\t3\t"+str(singlet_substitutions)+tab+str(singlet_substitutions1)+tab+str(singlet_substitutions2)+newline)
                f.write("Double_substitutions"+tab+run_id+"\t3\t"+str(doublet_substitutions)+tab+str(doublet_substitutions1)+tab+str(doublet_substitutions2)+newline)
                f.write("Triple_substitutions"+tab+run_id+"\t3\t"+str(triplet_substitutions)+tab+str(triplet_substitutions1)+tab+str(triplet_substitutions2)+newline)
                f.write("Quadruple_substitutions"+tab+run_id+"\t3\t"+str(quadruplet_substitutions)+tab+str(quadruplet_substitutions1)+tab+str(quadruplet_substitutions2)+newline)
                f.write("Quintuple_substitutions"+tab+run_id+"\t3\t"+str(quintuplet_substitutions)+tab+str(quintuplet_substitutions1)+tab+str(quintuplet_substitutions2)+newline)
                f.write("Sextuple_substitutions"+tab+run_id+"\t3\t"+str(sextuplet_substitutions)+tab+str(sextuplet_substitutions1)+tab+str(sextuplet_substitutions2)+newline)
                f.write("Septuple_substitutions"+tab+run_id+"\t3\t"+str(septuplet_substitutions)+tab+str(septuplet_substitutions1)+tab+str(septuplet_substitutions2)+newline)
                f.write("Octuple_substitutions"+tab+run_id+"\t3\t"+str(octuplet_substitutions)+tab+str(octuplet_substitutions1)+tab+str(octuplet_substitutions2)+newline)
                f.write("Longer_substitutions"+tab+run_id+"\t3\t"+str(longer_substitutions)+tab+str(longer_substitutions1)+tab+str(longer_substitutions2)+newline)

                f.write("Insertion_events"+tab+run_id+"\t3\t"+str(any_insertions)+tab+str(any_insertions1)+tab+str(any_insertions2)+newline)		
                f.write("Single_insertion"+tab+run_id+"\t3\t"+str(singlet_insertions)+tab+str(singlet_insertions1)+tab+str(singlet_insertions2)+newline)		
                f.write("Double_insertion"+tab+run_id+"\t3\t"+str(doublet_insertions)+tab+str(doublet_insertions1)+tab+str(doublet_insertions2)+newline)		
                f.write("Triple_insertion"+tab+run_id+"\t3\t"+str(triplet_insertions)+tab+str(triplet_insertions1)+tab+str(triplet_insertions2)+newline)		
                f.write("Quadruple_insertions"+tab+run_id+"\t3\t"+str(quadruplet_insertions)+tab+str(quadruplet_insertions1)+tab+str(quadruplet_insertions2)+newline)
                f.write("Quintuple_insertions"+tab+run_id+"\t3\t"+str(quintuplet_insertions)+tab+str(quintuplet_insertions1)+tab+str(quintuplet_insertions2)+newline)
                f.write("Sextuple_insertions"+tab+run_id+"\t3\t"+str(sextuplet_insertions)+tab+str(sextuplet_insertions1)+tab+str(sextuplet_insertions2)+newline)
                f.write("Septuple_insertions"+tab+run_id+"\t3\t"+str(septuplet_insertions)+tab+str(septuplet_insertions1)+tab+str(septuplet_insertions2)+newline)
                f.write("Octuple_insertions"+tab+run_id+"\t3\t"+str(octuplet_insertions)+tab+str(octuplet_insertions1)+tab+str(octuplet_insertions2)+newline)
                f.write("Insertion_up_to_40"+tab+run_id+"\t3\t"+str(insertions_40)+tab+str(insertions_40_1)+tab+str(insertions_40_2)+newline)		
                f.write("Insertion_up_to_1000"+tab+run_id+"\t3\t"+str(insertions_1000)+tab+str(insertions_1000_1)+tab+str(insertions_1000_2)+newline)		
                f.write("Insertion_long"+tab+run_id+"\t3\t"+str(insertions_long)+tab+str(insertions_long1)+tab+str(insertions_long2)+newline)		

                f.write("Deletion_events"+tab+run_id+"\t3\t"+str(any_deletions)+tab+str(any_deletions1)+tab+str(any_deletions2)+newline)
                f.write("Single_deletion"+tab+run_id+"\t3\t"+str(singlet_deletions)+tab+str(singlet_deletions1)+tab+str(singlet_deletions2)+newline)		
                f.write("Double_deletion"+tab+run_id+"\t3\t"+str(doublet_deletions)+tab+str(doublet_deletions1)+tab+str(doublet_deletions2)+newline)		
                f.write("Triple_deletion"+tab+run_id+"\t3\t"+str(triplet_deletions)+tab+str(triplet_deletions1)+tab+str(triplet_deletions2)+newline)		
                f.write("Quadruple_deletions"+tab+run_id+"\t3\t"+str(quadruplet_deletions)+tab+str(quadruplet_deletions1)+tab+str(quadruplet_deletions2)+newline)
                f.write("Quintuple_deletions"+tab+run_id+"\t3\t"+str(quintuplet_deletions)+tab+str(quintuplet_deletions1)+tab+str(quintuplet_deletions2)+newline)
                f.write("Sextuple_deletions"+tab+run_id+"\t3\t"+str(sextuplet_deletions)+tab+str(sextuplet_deletions1)+tab+str(sextuplet_deletions2)+newline)
                f.write("Septuple_deletions"+tab+run_id+"\t3\t"+str(septuplet_deletions)+tab+str(septuplet_deletions1)+tab+str(septuplet_deletions2)+newline)
                f.write("Octuple_deletions"+tab+run_id+"\t3\t"+str(octuplet_deletions)+tab+str(octuplet_deletions1)+tab+str(octuplet_deletions2)+newline)
                f.write("Deletion_up_to_40"+tab+run_id+"\t3\t"+str(deletions_40)+tab+str(deletions_40_1)+tab+str(deletions_40_2)+newline)		
                f.write("Deletion_up_to_1000"+tab+run_id+"\t3\t"+str(deletions_1000)+tab+str(deletions_1000_1)+tab+str(deletions_1000_2)+newline)		
                f.write("Deletion_long"+tab+run_id+"\t3\t"+str(deletions_long)+tab+str(deletions_long1)+tab+str(deletions_long2)+newline)		

                f.write("Single_insertion_in_polymers"+tab+run_id+"\t3\t"+str(singlet_s_insertions)+tab+str(singlet_s_insertions1)+tab+str(singlet_s_insertions2)+newline)
                f.write("Double_insertion_in_polymers"+tab+run_id+"\t3\t"+str(doublet_s_insertions)+tab+str(doublet_s_insertions1)+tab+str(doublet_s_insertions2)+newline)
                f.write("Triple_insertion_in_polymers"+tab+run_id+"\t3\t"+str(triplet_s_insertions)+tab+str(triplet_s_insertions1)+tab+str(triplet_s_insertions2)+newline)
                f.write("Quadruple_insertions_in_polymers"+tab+run_id+"\t3\t"+str(quadruplet_s_insertions)+tab+str(quadruplet_s_insertions1)+tab+str(quadruplet_s_insertions2)+newline)
                f.write("Quintuple_insertions_in_polymers"+tab+run_id+"\t3\t"+str(quintuplet_s_insertions)+tab+str(quintuplet_s_insertions1)+tab+str(quintuplet_s_insertions2)+newline)
                f.write("Sextuple_insertions_in_polymers"+tab+run_id+"\t3\t"+str(sextuplet_s_insertions)+tab+str(sextuplet_s_insertions1)+tab+str(sextuplet_s_insertions2)+newline)
                f.write("Septuple_insertions_in_polymers"+tab+run_id+"\t3\t"+str(septuplet_s_insertions)+tab+str(septuplet_s_insertions1)+tab+str(septuplet_s_insertions2)+newline)
                f.write("Octuple_insertions_in_polymers"+tab+run_id+"\t3\t"+str(octuplet_s_insertions)+tab+str(octuplet_s_insertions1)+tab+str(octuplet_s_insertions2)+newline)
                f.write("Longer_insertion_in_polymers"+tab+run_id+"\t3\t"+str(longer_s_insertions)+tab+str(longer_s_insertions1)+tab+str(longer_s_deletions2)+newline)
                f.write("Single_deletion_in_polymers"+tab+run_id+"\t3\t"+str(singlet_s_deletions)+tab+str(singlet_s_deletions1)+tab+str(singlet_s_deletions2)+newline)
                f.write("Double_deletion_in_polymers"+tab+run_id+"\t3\t"+str(doublet_s_deletions)+tab+str(doublet_s_deletions1)+tab+str(doublet_s_deletions2)+newline)
                f.write("Triple_deletion_in_polymers"+tab+run_id+"\t3\t"+str(triplet_s_deletions)+tab+str(triplet_s_deletions1)+tab+str(triplet_s_deletions2)+newline)
                f.write("Quadruple_deletions_in_polymers"+tab+run_id+"\t3\t"+str(quadruplet_s_deletions)+tab+str(quadruplet_s_deletions1)+tab+str(quadruplet_s_deletions2)+newline)
                f.write("Quintuple_deletions_in_polymers"+tab+run_id+"\t3\t"+str(quintuplet_s_deletions)+tab+str(quintuplet_s_deletions1)+tab+str(quintuplet_s_deletions2)+newline)
                f.write("Sextuple_deletions_in_polymers"+tab+run_id+"\t3\t"+str(sextuplet_s_deletions)+tab+str(sextuplet_s_deletions1)+tab+str(sextuplet_s_deletions2)+newline)
                f.write("Septuple_deletions_in_polymers"+tab+run_id+"\t3\t"+str(septuplet_s_deletions)+tab+str(septuplet_s_deletions1)+tab+str(septuplet_s_deletions2)+newline)
                f.write("Octuple_deletions_in_polymers"+tab+run_id+"\t3\t"+str(octuplet_s_deletions)+tab+str(octuplet_s_deletions1)+tab+str(octuplet_s_deletions2)+newline)
                f.write("Longer_deletion_in_polymers"+tab+run_id+"\t3\t"+str(longer_s_deletions)+tab+str(longer_s_deletions1)+tab+str(longer_s_deletions2)+newline)

                if sliding_insertions > 0:
                        f.write("Insertion_events_in_polymers"+tab+run_id+"\t3\t"+str(sliding_insertions)+tab+str(sliding_insertions1)+tab+str(sliding_insertions2)+newline)		
                if sliding_deletions > 0:
                        f.write("Deletion_events_in_polymers"+tab+run_id+"\t3\t"+str(sliding_deletions)+tab+str(sliding_deletions1)+tab+str(sliding_insertions2)+newline)		
 
                if soft_clipped_events > 0:
                        f.write("Soft_clipped_events"+tab+run_id+"\t1\t"+str(soft_clipped_events)+newline)		
                        f.write("Soft_clipped_bases"+tab+run_id+"\t1\t"+str(soft_clipped_bases)+newline)		
                if hard_clipped_events > 0:
                        f.write("Hard_clipped_events"+tab+run_id+"\t1\t"+str(hard_clipped_events)+newline)		
                        f.write("Hard_clipped_bases"+tab+run_id+"\t1\t"+str(hard_clipped_bases)+newline)		

                if intron_N > 0:
                        f.write("Intron"+tab+run_id+"\t3\t"+str(intron_N)+tab+str(intron_N_1)+tab+str(intron_N_2)+newline)		
                if intron_n > 0:
                        f.write("Intron_n"+tab+run_id+"\t3\t"+str(intron_n)+tab+str(intron_n_1)+tab+str(intron_n_2)+newline)		
                        
                if prefix_S > 0:
                        f.write("Unaligned_prefix_1_16"+tab+run_id+"\t3\t"+str(prefix_S)+tab+str(prefix_S_1)+tab+str(prefix_S_2)+newline)		
                if prefix_M > 0:
                        f.write("Unaligned_prefix_17_100"+tab+run_id+"\t3\t"+str(prefix_M)+tab+str(prefix_M_1)+tab+str(prefix_M_2)+newline)		
                if prefix_L > 0:
                        f.write("Unaligned_prefix_long"+tab+run_id+"\t3\t"+str(prefix_L)+tab+str(prefix_L_1)+tab+str(prefix_L_2)+newline)		

                if suffix_S > 0:
                        f.write("Unaligned_suffix_1_16"+tab+run_id+"\t3\t"+str(suffix_S)+tab+str(suffix_S_1)+tab+str(suffix_S_2)+newline)		
                if suffix_M > 0:
                        f.write("Unaligned_suffix_17_100"+tab+run_id+"\t3\t"+str(suffix_M)+tab+str(suffix_M_1)+tab+str(suffix_M_2)+newline)		
                if suffix_L > 0:
                        f.write("Unaligned_suffix_long"+tab+run_id+"\t3\t"+str(suffix_L)+tab+str(suffix_L_1)+tab+str(suffix_L_2)+newline)		


        if every_base_flag:
                keys = ["AA","TT","GG","CC","AG","TC","GA","CT","AT","TA","GC","CG","AC","TG","GT","CA"]  
        else:
                keys = ["AG","TC","GA","CT","AT","TA","GC","CG","AC","TG","GT","CA"]

        for k in keys:
                if not k in base_matrix_total : # jfm 9/28/2012 some types may not occur
                        base_matrix_total[k] = 0
                if not k in base_matrix_total1 : # jfm 9/28/2012 some types may not occur
                        base_matrix_total1[k] = 0
                if not k in base_matrix_total2 : # jfm 9/28/2012 some types may not occur
                        base_matrix_total2[k] = 0
                if view_type  == "tsv":
                        f.write("Sub:" + k +tab+run_id+ "\t3\t" + str(base_matrix_total[k]) + tab + str(base_matrix_total1[k]) + tab + str(base_matrix_total2[k])+newline)

###############################################################################################	
## Export the multi alignments

def exportMultiAlignments(f,run_id,view_type,aliHisto,aliHisto1,aliHisto2) :

	nmax = 0
	for k in aliHisto.keys() :
		if k > nmax and aliHisto[k] + aliHisto1[k] + aliHisto2[k] > 0 :
			nmax = k
	for k in range (nmax) :
		if view_type  == "tsv":
			f.write(("MultiAli:%03d" % (k+1)) +tab+run_id+ "\t3")
			f.write(tab + str( aliHisto[k+1] + aliHisto1[k+1] + aliHisto2[k+1] ))
			f.write(tab + str( aliHisto1[k+1]))
			f.write(tab + str( aliHisto2[k+1]) +newline )

##############################################################################################		
			
def registerAli (ali,ali1,ali2,aliHisto,aliHisto1,aliHisto2) :

	n = ali 
	if n > 100 :
		n = 100
	aliHisto[n] += 1
		
	n = ali1
	if n > 100 :
		n = 100
	aliHisto1[n] += 1
			
	n = ali2
	if n > 100 :
		n = 100
	aliHisto2[n] += 1

###############################################################################################		
###############################################################################################	
#
# Main program
#
###############################################################################################	
# get command line options and arguments

phase = 0
nphase = 0
nReadsGiven = 0
nBasesGiven = 0
sam_file = ""
fasta_file = ""
out_prefix = ""
input_file = ""
file_type = "SAM"
doMerge = 0
minAli = 0 
view_type = "tsv"
group_file = ""
run_id = ""
multi_file = 0
every_base_flag = False # jfm 10/12/2012 When False, skip base by base comparison when read == reference (Don't count AA,TT,CC,GG)

try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], '?ho:i:I:f:h:r:ameBSg:',['BAM','SAM','SAMGZ','SAMSORTED','SAMSORTEDGZ', 'magic','merge', 'split', 'view=','group=','minAli=','fasta=','run=','nreads=','nbases=', 'help'])
except getopt.GetoptError as err:
        print (str(err))
        print ("try --help")
        sys.exit(2)
        
if len(args) > 0:
        print ("unknown argument "+args[0])
        print ("try --help")
        sys.exit(2)

for o, a in opts:
        if o  == "-?" or o == "-h" or o == "--help":
                usage()
        elif o == "-o":
                out_prefix = a
        elif o == "-i":
                input_file = a
        elif o == "-r" or o == "--run":
                run_id = a
        elif o == "-m" or o == "--merge":
                doMerge = 1
                phase = 2
                nphase = nphase + 1 
        elif o == "-f" or  o == "--fasta":
                fasta_file = a
        elif o == "-a" or o == "--magic":
                file_type = "Magic"
                phase = 1
                nphase = nphase + 1 
        elif o == "-S" or o == "--SAM":
                file_type = "SAM"
                phase = 1
                nphase = nphase + 1 
        elif o == "--SAMGZ":
                file_type = "SAMGZ"
                phase = 1
                nphase = nphase + 1 
        elif o == "--SAMSORTED":
                file_type = "SAMSORTED"
                phase = 1
                nphase = nphase + 1 
        elif o == "--SAMSORTEDGZ":
                file_type = "SAMSORTEDGZ"
                phase = 1
                nphase = nphase + 1 
        elif o == "-B" or  o == "--BAM":
                file_type = "BAM"
                phase = 1
                nphase = nphase + 1 
        elif o == "--view":
                view_type = a
                phase = 3
                nphase = nphase + 1 
        elif o == "--nreads":
                nReadsGiven = int(a)
        elif o == "--nbases":
                nBasesGiven = int(a)
        elif o == "--minAli":
                minAli = int(a)
        elif o == "-g" or o == "--group":
                group_file = a
        elif o == "-e":
                every_base_flag = True
        elif o == "--split":
                multi_file = 1


###############################################################################################	
# validate the parameters 

if nphase == 0:
	print ("# no action specified, please try --help")
	sys.exit(1)
if nphase > 1:
	print ("# incompatible arguments, only one phase should be specified, please try --help")
	sys.exit(2)
	
if phase == 1:
	if run_id == "":
		print ("# In phase 1, it is mandatory to specify --run run_id, please try --help")
		sys.exit(3)
	if fasta_file == "":
		print ("# In phase 1, it is mandatory to specify --fasta target.fa, please try --help")
		sys.exit(4)
elif phase == 2:
	if group_file == "" and run_id == "":
		print ("# In phase 2, it is mandatory to specify either --run or --group,  please try --help")
		sys.exit(5)
elif phase == 3:
        if run_id == "":
                run_id = "any"
        if 0 and out_prefix == "":
                print ("# In phase 3, it is mandatory to specify -o file_prefix, please try --help")
                sys.exit(6)
        if view_type != "table"	and view_type != "html" :
                print ("# In phase 3, currently the only known formats are 'table' and 'html', please try --help")
                print ("# for the moment only the phase2 output is completed")
                print ("# the table output is very preliminary and incomplete")
                print ("# if you sort it, it will look a bit nicer")
                print ("# we will have better displays in the next few days")
                sys.exit(6)

###############################################################################################
# Select the input stream
gg =  """ gawk -F '\\t' '{gsub("*","0",$7);printf("%s",$1);for(i=2;i<=NF;i++)printf("\\t%s",$i);printf("\\n");}' """

if input_file == "":
        input_stream = os.popen( " sort -k 1,1 ")
        # sys.stdin
else:
        if file_type == "BAM":
                input_stream = os.popen( "samtools view -h " + input_file + " | " + gg + " | sort -T . -k 1,1 ")
        elif file_type == "SAMSORTED":
                input_stream = os.popen( "cat " + input_file + " | " + gg )
        elif file_type == "SAMSORTEDGZ":
                input_stream = os.popen( "gunzip -c " + input_file +  " | " + gg )
        elif file_type == "SAMGZ":
                input_stream = os.popen( "gunzip -c " + input_file +  " | " + gg + " | sort -T . -k 1,1 ")
        else:
                input_stream = os.popen( "cat " + input_file )

###############################################################################################		
# Define paths and filenames

if out_prefix == "" and phase < 3:
        tsv_file = ""
else:
        tsv_file = out_prefix + ".aliqc.tsv"
summary_file = out_prefix + ".summary.tsv"
any_event_file = out_prefix + ".mismatch_rate_hist.tsv"
mismatch_type_file = out_prefix + ".mismatch_type_hist.tsv"
aligned_length_file = out_prefix + ".match_length_hist.tsv"
multiple_alignment_file = out_prefix + ".multiple_alignment_hist.tsv"
mismatch_counts_per_cycle_file = out_prefix + ".mismatch_counts_per_cycle_file.tsv"
mismatch_per_cycle_per_kb_aligned_file = out_prefix + ".mismatch_per_cycle_per_kb_aligned.tsv"
mismatch_histo_and_types_file = out_prefix + ".mismatch_histo_and_types.tsv"
aligned_reads_multiplicity_length_aligned_file = out_prefix + ".aligned_reads_multiplicity_length_aligned.tsv"

indel_type_file = out_prefix + ".indel_type.tsv"
err_pos_file = out_prefix + ".err_pos.tsv"

tab = "\t"
newline = "\n"
 
# mieg: i change the interface, the user is now required to provide a prefix which
# may include the name of an existing
# directory, he may also choose to export to stdout
# we could be more fancy, split the prefix and create the hiearchy
# for example if the prefix is a/b/c we would mkdir a a/b and export in a/b/c.aliqc.tsv
# but this may be invasive if we happen to create subdir in simlinked directories
# On a parallelized architecture, i am afraid of creating havock on people's disk
# i'd rather let them take their responsabilities outside of our code !
# Create the output directory
# if not os.path.exists(out_prefix):
#	os.mkdir(out_prefix)
	
###############################################################################################		
# in merge case, parse the group file, cumulate and reexport
if phase == 2:
	phase2(run_id,input_stream, out_prefix,nReadsGiven,nBasesGiven)
	sys.exit(0)

###############################################################################################		
# in merge case, parse the group file, cumulate and reexport
if phase == 3:
	phase3(input_stream, view_type, out_prefix,nReadsGiven,nBasesGiven,multi_file)
	sys.exit(0)

###############################################################################################		
# In phase 1, we want to process a BAM file, using the HTseq library
# by having import here in the code, we can run phase2, or phase3 even if HTseq is not available
import HTSeq
print ("# HTSeq succesfully imported")
###############################################################################################		
# Define histograms as python dictionaries

alignments_per_read = {}
alignments_per_read1 = {}
alignments_per_read2 = {}

match_length_hist = {}
match_length_hist1 = {}
match_length_hist2 = {}

any_event_hist = {}
any_event_hist1 = {}
any_event_hist2 = {}

base_matrix_total = {}
base_matrix_total1 = {}
base_matrix_total2 = {}

valid_bases = "ATGC"     # mieg reorder ATGC
for i in valid_bases:
	for j in valid_bases:	# Optimize base matrix totals by creating keys and initializing (jfm 10/12/2012 opt2)
		c = i + j
		base_matrix_total[c] = 0
		base_matrix_total1[c] = 0
		base_matrix_total2[c] = 0

delete_type_hist = {}
delete_type_hist1 = {}
delete_type_hist2 = {}

insert_type_hist = {}
insert_type_hist1 = {}
insert_type_hist2 = {}

err_pos = {}
err_pos1 = {}
err_pos2 = {}
err_N1 = {}
err_N2 = {}

###############################################################################################		
# read the entire reference genome file into a dictionary
print ("# Reading fasta sequences from "+fasta_file+"... "+ time.asctime(time.localtime(time.time())))
seq = dict( (s.name, s) for s in HTSeq.FastaReader( fasta_file ))
print ("# Done. " + time.asctime(time.localtime(time.time())))
print ("# Processing file "+input_file)
print ("# Output prefix is "+out_prefix)
print ()

###############################################################################################
# Process a Magic file
# in this format, the big alignments files have already been preprocessed in C
# we are simply parsing an ace file for class Ali
# as usual, we ignore the run names and add up corresponding values, so the code is simple and fast

# parseExportMagicFile (run_id,input_stream, out_prefix)
# exit (0)

###############################################################################################
# Process the BAM/SAM file
# The BAM file countains records
record_count = 0
# Most correspond to alignments
old = ""
reads_in_file = 0
reads_1_in_file = 0
reads_2_in_file = 0
alignment_count = 0
# number of aligned bases
read_bp = 0
read_1_bp = 0
read_2_bp = 0
readU_bp = 0
readU_1_bp = 0
readU_2_bp = 0
ali_bp = 0
ali_1_bp = 0
ali_2_bp = 0
match_base = 0
sub_base = 0
del_base = 0
ins_base = 0
match_base_1 = 0
sub_base_1 = 0
del_base_1 = 0
ins_base_1 = 0
match_base_2 = 0
sub_base_2 = 0
del_base_2 = 0
ins_base_2 = 0

# in paired-end protocols, we can decompose
read1_alignment_count = 0
read2_alignment_count = 0

# the remainder of the BAM file is a way to give the unaligned reads
unaligned_count = 0
read1_unaligned_count = 0
read2_unaligned_count = 0

# From here on, each read will be counted only once
reads_mapped_count = 0
read1_mapped_count = 0
read2_mapped_count = 0
# in stranded protocol, the stranding is usually opposite for the 2 reads and must be counted separatelly
reads_plus_strand_count = 0
reads_minus_strand_count = 0
read1_plus_strand_count = 0
read1_minus_strand_count = 0
read2_plus_strand_count = 0
read2_minus_strand_count = 0

# Each read contributes once to the mismatches and the pair statistics
perfect = 0
perfect1 = 0
perfect2 = 0

perfect_clipped = 0
perfect_clipped1 = 0
perfect_clipped2 = 0

hard_clipped_events = 0
hard_clipped_bases = 0
soft_clipped_events = 0
soft_clipped_bases = 0

intron_N = 0
intron_N_1 = 0
intron_N_2 = 0

intron_n = 0
intron_n_1 = 0
intron_n_2 = 0

prefix_S = 0
prefix_S_1 = 0
prefix_S_2 = 0

suffix_S = 0
suffix_S_1 = 0
suffix_S_2 = 0

prefix_M = 0
prefix_M_1 = 0
prefix_M_2 = 0

suffix_M = 0
suffix_M_1 = 0
suffix_M_2 = 0

prefix_L = 0
prefix_L_1 = 0
prefix_L_2 = 0

suffix_L = 0
suffix_L_1 = 0
suffix_L_2 = 0

any_substitutions = 0
any_substitutions1 = 0
any_substitutions2 = 0

any_insertions = 0
any_insertions1 = 0
any_insertions2 = 0

any_deletions = 0
any_deletions1 = 0
any_deletions2 = 0

insertions_40 = 0
insertions_40_1 = 0
insertions_40_2 = 0

insertions_1000 = 0
insertions_1000_1 = 0
insertions_1000_2 = 0

insertions_long = 0
insertions_long1 = 0
insertions_long2 = 0

deletions_40 = 0
deletions_40_1 = 0
deletions_40_2 = 0

deletions_1000 = 0
deletions_1000_1 = 0
deletions_1000_2 = 0

deletions_long = 0
deletions_long1 = 0
deletions_long2 = 0

sliding_insertions = 0
sliding_insertions1 = 0
sliding_insertions2 = 0

sliding_deletions = 0
sliding_deletions1 = 0
sliding_deletions2 = 0

singlet_substitutions = 0
singlet_substitutions1 = 0
singlet_substitutions2 = 0

singlet_insertions = 0
singlet_insertions1 = 0
singlet_insertions2 = 0

singlet_deletions = 0
singlet_deletions1 = 0
singlet_deletions2 = 0

doublet_substitutions = 0
doublet_substitutions1 = 0
doublet_substitutions2 = 0

doublet_insertions = 0
doublet_insertions1 = 0
doublet_insertions2 = 0

doublet_deletions = 0
doublet_deletions1 = 0
doublet_deletions2 = 0

triplet_substitutions = 0
triplet_substitutions1 = 0
triplet_substitutions2 = 0

longer_substitutions = 0
longer_substitutions1 = 0
longer_substitutions2 = 0

triplet_insertions = 0
triplet_insertions1 = 0
triplet_insertions2 = 0

triplet_deletions = 0
triplet_deletions1 = 0
triplet_deletions2 = 0

quadruplet_substitutions = 0
quadruplet_substitutions1 = 0
quadruplet_substitutions2 = 0

quadruplet_insertions = 0
quadruplet_insertions1 = 0
quadruplet_insertions2 = 0

quadruplet_deletions = 0
quadruplet_deletions1 = 0
quadruplet_deletions2 = 0

quintuplet_substitutions = 0
quintuplet_substitutions1 = 0
quintuplet_substitutions2 = 0

quintuplet_insertions = 0
quintuplet_insertions1 = 0
quintuplet_insertions2 = 0

quintuplet_deletions = 0
quintuplet_deletions1 = 0
quintuplet_deletions2 = 0

sextuplet_substitutions = 0
sextuplet_substitutions1 = 0
sextuplet_substitutions2 = 0

sextuplet_insertions = 0
sextuplet_insertions1 = 0
sextuplet_insertions2 = 0

sextuplet_deletions = 0
sextuplet_deletions1 = 0
sextuplet_deletions2 = 0

septuplet_substitutions = 0
septuplet_substitutions1 = 0
septuplet_substitutions2 = 0

septuplet_insertions = 0
septuplet_insertions1 = 0
septuplet_insertions2 = 0

septuplet_deletions = 0
septuplet_deletions1 = 0
septuplet_deletions2 = 0

octuplet_substitutions = 0
octuplet_substitutions1 = 0
octuplet_substitutions2 = 0

octuplet_insertions = 0
octuplet_insertions1 = 0
octuplet_insertions2 = 0

octuplet_deletions = 0
octuplet_deletions1 = 0
octuplet_deletions2 = 0

singlet_s_insertions = 0
singlet_s_insertions1 = 0
singlet_s_insertions2 = 0

singlet_s_deletions = 0
singlet_s_deletions1 = 0
singlet_s_deletions2 = 0

doublet_s_insertions = 0
doublet_s_insertions1 = 0
doublet_s_insertions2 = 0

doublet_s_deletions = 0
doublet_s_deletions1 = 0
doublet_s_deletions2 = 0

triplet_s_insertions = 0
triplet_s_insertions1 = 0
triplet_s_insertions2 = 0

triplet_s_deletions = 0
triplet_s_deletions1 = 0
triplet_s_deletions2 = 0

quadruplet_s_insertions = 0
quadruplet_s_insertions1 = 0
quadruplet_s_insertions2 = 0

quadruplet_s_deletions = 0
quadruplet_s_deletions1 = 0
quadruplet_s_deletions2 = 0

quintuplet_s_insertions = 0
quintuplet_s_insertions1 = 0
quintuplet_s_insertions2 = 0

quintuplet_s_deletions = 0
quintuplet_s_deletions1 = 0
quintuplet_s_deletions2 = 0

sextuplet_s_insertions = 0
sextuplet_s_insertions1 = 0
sextuplet_s_insertions2 = 0

sextuplet_s_deletions = 0
sextuplet_s_deletions1 = 0
sextuplet_s_deletions2 = 0

septuplet_s_insertions = 0
septuplet_s_insertions1 = 0
septuplet_s_insertions2 = 0

septuplet_s_deletions = 0
septuplet_s_deletions1 = 0
septuplet_s_deletions2 = 0

octuplet_s_insertions = 0
octuplet_s_insertions1 = 0
octuplet_s_insertions2 = 0

octuplet_s_deletions = 0
octuplet_s_deletions1 = 0
octuplet_s_deletions2 = 0

longer_s_insertions = 0
longer_s_insertions1 = 0
longer_s_insertions2 = 0

longer_s_deletions = 0
longer_s_deletions1 = 0
longer_s_deletions2 = 0

bad_base_count = 0
proper_pair_count = 0
pair_mapped_count = 0

ali = 0
ali1 = 0
ali2 = 0

primary = 0
primary1 = 0
primary2 = 0

aliHisto = {}
aliHisto1 = {}
aliHisto2 = {}

for n in range(101):
	aliHisto[n] = 0
	aliHisto1[n] = 0
	aliHisto2[n] = 0

# Process each record in the SAM file
print ("BAM/SAM processing starts " + time.asctime(time.localtime(time.time())))

for a in HTSeq.SAM_Reader( input_stream ):

        record_count += 1
        if  0 and 	record_count > 10000:
                break
        
        if minAli > 0:
                bpAli = 0
                if a.cigar:
                        for i in range(0, len(a.cigar)):
                                type = a.cigar[i].type
                                if type == "M" or type == "=" or type == "X": # Alignment match (can be a sequence match or mismatch)
                                        bpAli += a.cigar[i].size
                if bpAli < minAli :
                        continue
                                       
                                       
        if a.aligned:
                alignment_count += 1

                # if a.not_primary_alignment: # jfm 10/18/2012 Skip secondary alignments
                # continue

                first_alignment = False	# jfm 10/1/2012 Only do statistics on the first alignment for each end

                # count read1 and read2
                if 1 :
                        if a.read.name != old :
                                old = a.read.name
			
                                registerAli (ali,ali1,ali2,aliHisto,aliHisto1,aliHisto2)
				
                                ali = 0
                                ali1 = 0
                                ali2 = 0

                                primary = 0
                                primary1 = 0
                                primary2 = 0

                        if a.paired_end:
                                if a.pe_which == "first":
                                        read1_alignment_count += 1
                                        ali1 += 1
                                        if not a.not_primary_alignment:
                                                primary1 += 1
                                                if primary1 == 1:
                                                        first_alignment = True # jfm 10/1/2012 Only do statistics on the first alignment for each end
                                elif a.pe_which == "second":
                                        read2_alignment_count += 1
                                        ali2 += 1
                                        if not a.not_primary_alignment:
                                                primary2 += 1
                                                if primary2 == 1:
                                                        first_alignment = True # jfm 10/1/2012 Only do statistics on the first alignment for each end
                        else:
                                ali += 1
                                if not a.not_primary_alignment:
                                        primary += 1
                                        if primary == 1:
                                                first_alignment = True

                # mieg: do basic alignment statistics ONLY for first alignment
                # we want to count the aligned reads and their mismatches not the alignments,
                # otherwise a bugged code aligning each read at 10000 different places will
                # seem to align more sequences than a correct aligner
		# jfm 10/1/2012 Do statistics on the first alignment for each end
                if ( first_alignment ):
                        if calc_seqc(a) == 0:
                                continue
                        reads_in_file += 1
                        reads_mapped_count += 1
                        # get paired end counts
                        if a.paired_end:
		
                                # first and second fragment alignments  
                                if a.pe_which == "first":
                                        reads_1_in_file += 1		    
                                        read1_mapped_count += 1
                                elif a.pe_which == "second":
                                        reads_2_in_file += 1		    
                                        read2_mapped_count += 1
				
                                # with itself and mate mapped
                                if a.mate_aligned:
                                        pair_mapped_count += 1

                                # proper pairs
                                if a.proper_pair:
                                        proper_pair_count += 1
                

        else:
                unaligned_count += 1
                reads_in_file += 1
                k1 = len(a.read)
                if k1 == 1 :
                        k1 = 0  # it is not interesting to count the length of a pseudo-read == "*"
                readU_bp += k1

                if a.paired_end:
                        if a.pe_which == "first":
                                read1_unaligned_count += 1
                                reads_1_in_file += 1
                                readU_1_bp += k1
                        elif a.pe_which == "second":
                                read2_unaligned_count += 1
                                reads_2_in_file += 1
                                readU_2_bp += k1
	# del(a)

registerAli (ali,ali1,ali2,aliHisto,aliHisto1,aliHisto2)

# Debug
# print ("# Bad bases:", bad_base_count)

# Print (summary file
if 1 :
        if tsv_file == "":
                f = sys.stdout
        else:
                f = open(tsv_file, "w")

        exportSummary(f,run_id,view_type,nReadsGiven,nBasesGiven)
        exportAlignedLngths(f,run_id,view_type,match_length_hist)
        exportAnyEvent(f,run_id,view_type,any_event_hist)
        exportSubstitutions(f,run_id,view_type,base_matrix_total,base_matrix_total1,base_matrix_total2)
        exportDeletions(f,run_id,view_type,delete_type_hist, delete_type_hist1, delete_type_hist2)
        exportInsertion(f,run_id,view_type,insert_type_hist,insert_type_hist1,insert_type_hist2)
        exportMultiAlignments(f,run_id,view_type,aliHisto,aliHisto1,aliHisto2)
        exportMismatchPositions(f,run_id,view_type,err_pos,err_pos1,err_pos2,err_N1,err_N2)

print ("processing done " + time.asctime(time.localtime(time.time())))

########################################################################################
########################################################################################
# Process a Magic file
# in this format, the big alignments files have already been preporcessed in C
# we are simply parsing an ace file for class Ali
# as usual, we ignore the run names and add up corresponding values, so the code is simple and fast


def parseExportMagicFile (run_id,input_stream, out_prefix) :

        exit (0)

########################################################################################
########################################################################################
########################################################################################
	
