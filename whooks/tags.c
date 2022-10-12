/*  File: tags.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: Kernel file
 *              Edit it if you add a new tag, or do the lexaddkey in your own code.
 *              
 * Exported functions:
 *              tagInit(), called from lexsubs.c
 * HISTORY:
 * Last edited: Jul 16 18:26 2001 (edgrif)
 * Created: Mon Aug 29 17:09:21 1994 (mieg)
 *-------------------------------------------------------------------
 */


/* $Id: tags.c,v 1.10 2012/03/30 15:18:49 mieg Exp $  */


/* VERY IMPORTANT,
  Most tags are now define as 
     extern KEY
  in the file whooks/tags.h.
  These must be allocated here and initialised dynamically.
  Therefore any addition in tags.wrm must be matched by a
  declaration and an initialisation here.

  Note that this is not necessary;
  a new application can very well define its own tags statically
  in its own file by doing a lexaddkey there, without ever touching 
  the present file  or the tags.h file.
  
  Note also that the systags follow a different system, they are
  #defined to numbers in systags.wrm and strictly initialised
  this is because they are needed at bootstrap before reading
  in the System lexique.
  */

#include "acedb.h"
#include "lex.h"
#include "dbpath.h"
#include "whooks/sysclass.h"
#include "whooks/tags.h"

KEY
  _1_does_not_include_2,
  _1_includes_2,
  _2_point,
  _A1_labelling,
  _A_non_B,
  _Abstract,
  _Address,
  _Aligned,
  _Aligned_into,
  _Allele,
  _Allele_designation,
  _Alu_segment,
  _Amber,
  _Approximate_Match_to,
  _Assembly_repeat,
  _Assembly_tags,
  _Ace_mbly_tags,
  _Author,
  _Autopos,
  _B_non_A,
  _Back_one,
  _Backcross,
  _Balances,
  _Bands,
  _Bitnet,
  _Brief_identification,
  _CDS,
  _CDS_predicted_by,
  _CGC,
  _Calc,
  _Calc_distance,
  _Calc_lower_conf,
  _Calc_upper_conf,
  _Calculation,
  _Calculation_1,
  _Calculation_2,
  _Canon_for_cosmid,
  _Canonical_for,
  _Chromosome,
  _Cleavage,
  _Clone,
  _Subclone,
  _Clone_as_locus,
  _Clone_inside,
  _Clone_outside,
  _Cloning_vector,
  _Coding,
  _Cold_sensitive,
  _Columns,
  _Combined,
  _Comment,
  _Company,
  _Complementation_data,
  _Complex_mixed,
  _Composite, 
  _Compound,
  _Compression,
  _Contained_in,
  _Contains,
  _Contents,
  _Contig,
  _Contig9,
  _Contiguous,
  _Correct_name,
  _Corresponding_DNA,
  _Corresponding_peptide,
  _Cosmid_grid,
  _Cosmid_vector,
  _Covers,
  _DB_remark,
  _DB_searched,
  _DNA,
  _DNA_homol,
  _Deletion,
  _Description,
  _Designating_laboratory,
  _Detection_method,
  _Df_Dup,
  _Direct,
  _Display,
  _Distance,
  _Does_not_hybridize_to,
  _Dom_let,
  _Dom_one,
  _Dom_selected,
  _Dom_semi,
  _Dominant,
  _Duplication,
  _E_mail,
  _End_not_found,
  _Error,
  _Error_scale,
  _Exact_Match_to,
  _Experiment,
  _Expression_construct,
  _Extent,
  _Fax,
  _Features,
  _FingerPrint,
  _Finished,
  _Flag,
  _Float_Params,
  _Free_dup,
  _Freezer,
  _From,
  _From_Author,
  _From_Laboratory,
  _From_Library,
  _From_left_end,
  _From_map,
  _Full,
  _Full_name,
  _Function,
  _Funny_Match_to,
  _Gel_Number,
  _Gel_length,
  _Gene,
  _GeneId,
  _Gene_class,
  _Gene_classes,
  _General,
  _General_remark,
  _Genetics,
  _Genomic_Canonical,
  _Genotype,
  _Gridded,
  _Heat_sensitive,
  _Homol,
  _Hybrid_cell_line,
  _Hybridizes_to,
  _In_Situ,
  _In_pool,
  _Inherits,
  _Inside,
  _Int_Params,
  _Internet,
  _Interval,
  _Intragenic_revertant_of_dominant,
  _Inverted,
  _Is_buried_under,
  _Isolation,
  _Isoschizomers,
  _Journal,
  _Keyword,
  _Laboratory,
#if !defined(MACINTOSH)
  _Layout,
#endif
  _Left,
  _Length,
  _Library,
  _Lines_at,
  _Linkage,
  _LiquidN2,
  _Location,
  _Locus,
  _LocusId,
  _Locus_1,
  _Locus_2,
  _Locus_A,
  _Locus_B,
  _Loop,
  _Mail,
  _Males,
  _Map,
  _Mapper,
  _Mapping_data,
  _Maps_with,
  _Marker_locus,
  _Matching_Genomic,
  _Matching_cDNA,
  _Maternal,
  _Max,
  _Min,
  _Min_score,
  _Minus70,
  _Molecular_information,
  _Multi_point,
  _Multiplet,
  _Muscle,
  _Mutagen,
  _N_gaps,
  _Name,
  _Negative_clone,
  _Negative_locus,
  _Negative_probe,
  _Nick_name,
  _No_overlap,
  _No_stagger,
  _Offset,
  _Old_CGC_distance,
  _Old_CGC_results,
  _Oligo,
  _One_all,
  _One_let,
  _One_recombinant,
  _Origin,
  _Other_name,
  _Outside,
  _Overhang,
  _Overlap,
  _PCR_product,
  _PCR_remark,
  _Page,
  _Paper,
  _Parameters,
  _Paternal,
  _Pattern,
  _Pep_homol,
  _Peptide,
  _Percent_Identity,
  _Phenotype,
  _Phone,
  _Polymorphism,
  _Position,
  _Positive_clone,
  _Positive_locus,
  _Positive_pool_probe,
  _Positive_probe,
  _Possible_exon,
  _Precursor,
  _Primers,
  _Procedure,
  _Processed_mRNA,
  _Product,
  _Properties,
  _Protocol,
  _Pseudogene,
  _Publisher,
  _Qualifier,
  _RFLP,
  _RNA,
  _Rearrangement,
  _Rearrangement_1,
  _Rearrangement_2,
  _Recessive,
  _Recs_all,
  _Reference,
  _Reference_Allele,
  _Reference_strain,
  _Refers_to,
  _Rejected,
  _Rejected_from,
  _Rejected_Reads,
  _Related_Sequence,
  _Remark,
  _Repeats,
  _Representative,
  _Results,
  _Right,
  _Row,
  _Score,
  _Search_Method,
  _Selected,
  _Selected_loci,
  _Selected_trans,
  _SELF,
  _Semi_dominant,
  _Seq_length,
  _Sequence,
  _Sequencing_vector,
  _Sex_cis,
  _Sex_full,
  _Sex_one,
  _Shotgun,
  _Significant_bases,
  _Similarity,
  _Site,
  _Source,
  _Source_Exons,
  _Space_at,
  _Species,
  _Staff,
  _Start_not_found,
#if !defined(MACINTOSH)
  _Status,
#endif
  _Stop,
  _Strain,
  _Strain_designation,
  _Strictly_Maternal,
  _Structure,
  _Subsequence,
  _TATA_signal,
  _TSL,
  _TSL_site,
  _Tandem,
  _Temperature,
  _Temperature_sensitive,
  _Tested,
  _Text_Params,
  _Title,
  _To_right_end,
  _Trace_quality,
  _Excellent_upto,
  _Good_upto,
  _Fair_upto,
  _Transcript,
  _Translocation,
  _Transposon,
  _Transposon_insertion,
  _Type,
  _Type_1,
  _Type_2,
  _Unit_Length,
  _Unprocessed_mRNA,
  _Variant,
  _Variant_of,
  _Vaxmap,
  _Vector,
  _Version,
  _Volume,
  _Weak,
  _Well_ordered,
  _With_Maternal_Effect,
  _Y_remark,
  _Year,
  _cDNA,
  _mRNA,
  _mat_peptide,
  _misc_feature,
  _misc_signal,
  _modified_base,
  _mutation,
  _old_sequence,
  _pMap,
  _polyA_signal,
  _polyA_site,
  _promoter,
  _rRNA,
  _repeat_region,
  _repeat_unit,
  _sig_peptide,
  _snRNA,
  _tRNA,
  _A_Repeat,
  _Note,
  _Subpool,
  _Colour,
  _ORF,
  _ATG,
  _Splice3,
  _Splice5,
  _Coding_seg,
  _Foreign_Reference,
  _Virtual_row,
  _In_grid,
  _Chrom_Band,
  _Centromere,
  _Dark,
  _NOR,
  _Inherits_from,
  _Fate,
  _Lineage,
  _Parent,
  _Daughter,
  _Group,
  _In_group,
  _Group_member,
  _Neurodata,
  _Send,
  _Send_joint,
  _Receive,
  _Receive_joint,
  _Gap_junction,
  _Contact,
  _Lineage_name,
  _Equivalence_fate,
  _Equivalence_origin,
  _Embryo_division_time,
  _Clone_left_end,
  _Clone_right_end,
  _Derived_from,
  _Derivative,
  _Isolated_for,
  _Wormpep,
  _Replaces,
  _Replaced_by,
  _Confidential_remark,
  _Submitted,
  _Overlap_right,
  _Overlap_left,
  _Repeat_consensus,
  _cDNA_EST,
  _Archived,
  _Lethal,
  _Hand_fixed,
  _Balancer,
  _Region,
  _Hand_verified,
  _YAC,
  _Cytogenetic,
  _Inside_YAC,
  _Inside_Fragment,
  _Bibliography,
  _Links,
  _Link,
  _Chimeric,
  _Non_Chimeric,
  _Size,
  _Interval_Mapping,
  _Main_Marker,
  _Ends,
  _Mapping,
  _p_Telomere,
  _q_Telomere,
  _Drawing,
  _Multi_Position,
  _Does_not_Contain,
  _Contained_Locus,
  _Exterior_Locus,
  _Inside_Rearr,
  _Outside_Rearr,
  _Overlaps,
  _Does_not_Overlap,
  _No_Overlap_Rearr,
  _Overlapping_Rearr,
  _Contained_Rearr,
  _Email,
  _Flipped,
  _Assembled_as_dummy,
  _Previous_contig,
  _Assembled_from,
  _Assembled_into,
  _Proposed,
  _BaseCall, 
  _Genetic,
  _Physical,
  _Restriction_Enzyme,
  _AA,
  _Clone_Grid,
  _Df_Dup_data,
  _Fragment,
  _Pages,
  _Medline_ID,
  _Lab_Location,
  _Match,
  _Method,
  _Motif,
  _MultiMap,
  _Multi_counts,
  _Multi_pt_data,
  _Pool,
  _Unpublished,
  _EC_enzyme,
  _Document,
  _Unit,
  _Enzyme,
  _GDB_id,
  _OMIM,
  _Probe,
  _ATCC_ref,
  _Add_date,
  _Annotation,
  _Approval_date,
  _Approved,
  _Assign_mode,
  _Availability,
  _CS,
  _Cosmid_clones,
  _Created,
  _Current_symbol,
  _Cytogenetic_location,
  _DNA_type,
  _Data_source,
  _EC_ref,
  _Editing_Date,
  _Excision_Data,
  _GDB,
  _GenBank_ref,
  _HGML_ref,
  _KW,
  _Links_to,
  _Max_heterozygosity,
  _Modification_date,
  _Modified,
  _OMIM_document,
  _Old_Num,
  _Originator,
  _Polymorphism_Data,
  _Previous_symbol,
  _Vector_type,
  _Primer,
  _SCF_File,
  _PCR,
  _Cor,
  _Neighbours,
  _Common_bands,
  _Band_Values,
  _Segment_Lengths,
  _Fragment_of,
  _ABI,
  _Wild_type,
  _ABI_Date,
  _ABI_Comment,
  _ABI_Machine,
  _Sample,
  _ABI_Analysis,
  _Run_Start,
  _Run_Stop,
  _Clipped,
  _Clipping,
  _Old_Clipping,
  _Hand_Clipping,
  _Vector_Clipping,
  _Symbol,
  _Anchor,
  _Band,
  _Gel,
  _Band_Lengths,
  _Otto,
  _Main_Locus,
  _STS,
  _Checked_individually,
  _Outside_YAC,
  _STS_out,
  _AluPCR,
  _by_AluPCR_Hybridisation,
  _by_Finger_Printing,
  _Kpn,
  _THE,
  _Somatic_hybrid,
  _Genethon,
  _Route_2,
  _Inside_YAC_Not_Verified,
  _Locus_Not_Verified,
  _Possibly_Contains_Locus,
  _Possibly_Inside_YAC,
  _Locus_Out,
  _Allelic_Variants,
  _Clinical_Synopses,
  _Cross_References,
  _Diagnosis,
  _Evolution,
  _Historical_Information,
  _Index_Terms,
  _Last_Edited,
  _Locus_Description,
  _MIM_Number,
  _Old_MIM_Number,
  _Population_Genetics,
  _References,
  _Treatment,
  _Main_Text,
  _Animal_Models,
  _map_location,
  _map_error,
  _LongText,
  _PCR_data,
  _Clone_data,
  _Oligo_DNA,
  _Allelep,
  _Allelem,
  _Gametep,
  _Gametem,
  _Tetrad,
  _Centromere_segregation ,
  _Genetic_code, 
  _Lethal_tested,

/* for Action class (added by P.Kocab) */

  _description,
  _used_in,
  _action_type,
  _internal,
  _ace_query,
  _ace_dump,
  _name_dump,
  _tace_command,
  _external,
  _compound,
  _action_element, 
  _params,
  _no_report_status,
  _read_ace,
  _synchro,
  _import_keyset,

/* for comparative maps (added by jld) */

  _Homology,
  _Doc,
  _Pairwise,
  _Fuzzy,
  _Symbol,
  _SMAP, 
  _S_Parent
;

void tagInit (void)
{

  lexaddkey("1_does_not_include_2", &_1_does_not_include_2, 0) ; 
  lexaddkey("1_includes_2", &_1_includes_2, 0) ; 
  lexaddkey("2_point", &_2_point, 0) ; 
  lexaddkey("A1_labelling", &_A1_labelling, 0) ; 
  lexaddkey("A_non_B", &_A_non_B, 0) ; 
  lexaddkey("Abstract", &_Abstract, 0) ; 
  lexaddkey("Address", &_Address, 0) ;
  lexaddkey("Aligned", &_Aligned, 0) ;
  lexaddkey("Aligned_into", &_Aligned_into, 0) ;
  lexaddkey("Allele", &_Allele, 0) ; 
  lexaddkey("Allele_designation", &_Allele_designation, 0) ; 
  lexaddkey("Alu_segment", &_Alu_segment, 0) ; 
  lexaddkey("Amber", &_Amber, 0) ; 
  lexaddkey("Approximate_Match_to", &_Approximate_Match_to, 0) ; 
  lexaddkey("Assembly_repeat", &_Assembly_repeat, 0) ; 
  lexaddkey("Assembly_tags", &_Assembly_tags, 0) ; 
  lexaddkey("Ace_mbly_tags", &_Ace_mbly_tags, 0) ; 
  lexaddkey("Author", &_Author, 0) ; 
  lexaddkey("Autopos", &_Autopos, 0) ; 
  lexaddkey("B_non_A", &_B_non_A, 0) ; 
  lexaddkey("Back_one", &_Back_one, 0) ; 
  lexaddkey("Backcross", &_Backcross, 0) ; 
  lexaddkey("Balances", &_Balances, 0) ; 
  lexaddkey("Bands", &_Bands, 0) ; 
  lexaddkey("Bitnet", &_Bitnet, 0) ; 
  lexaddkey("Brief_identification", &_Brief_identification, 0) ; 
  lexaddkey("CDS", &_CDS, 0) ; 
  lexaddkey("CDS_predicted_by", &_CDS_predicted_by, 0) ; 
  lexaddkey("CGC", &_CGC, 0) ; 
  lexaddkey("Calc", &_Calc, 0) ; 
  lexaddkey("Calc_distance", &_Calc_distance, 0) ; 
  lexaddkey("Calc_lower_conf", &_Calc_lower_conf, 0) ; 
  lexaddkey("Calc_upper_conf", &_Calc_upper_conf, 0) ; 
  lexaddkey("Calculation", &_Calculation, 0) ; 
  lexaddkey("Calculation_1", &_Calculation_1, 0) ; 
  lexaddkey("Calculation_2", &_Calculation_2, 0) ; 
  lexaddkey("Canon_for_cosmid", &_Canon_for_cosmid, 0) ; 
  lexaddkey("Canonical_for", &_Canonical_for, 0) ; 
  lexaddkey("Chromosome", &_Chromosome, 0) ; 
  lexaddkey("Cleavage", &_Cleavage, 0) ; 
  lexaddkey("Clone", &_Clone, 0) ; 
  lexaddkey("Subclone", &_Subclone, 0) ; 
  lexaddkey("Clone_as_locus", &_Clone_as_locus, 0) ; 
  lexaddkey("Clone_inside", &_Clone_inside, 0) ; 
  lexaddkey("Clone_outside", &_Clone_outside, 0) ; 
  lexaddkey("Cloning_vector", &_Cloning_vector, 0) ; 
  lexaddkey("Coding", &_Coding, 0) ; 
  lexaddkey("Cold_sensitive", &_Cold_sensitive, 0) ; 
  lexaddkey("Columns", &_Columns, 0) ; 
  lexaddkey("Combined", &_Combined, 0) ; 
  lexaddkey("Comment", &_Comment, 0) ; 
  lexaddkey("Company", &_Company, 0) ; 
  lexaddkey("Complementation_data", &_Complementation_data, 0) ; 
  lexaddkey("Complex_mixed", &_Complex_mixed, 0) ; 
  lexaddkey("Composite", &_Composite, 0) ; 
  lexaddkey("Compound", &_Compound, 0) ; 
  lexaddkey("Compression", &_Compression, 0) ; 
  lexaddkey("Contained_in", &_Contained_in, 0) ; 
  lexaddkey("Contains", &_Contains, 0) ; 
  lexaddkey("Contents", &_Contents, 0) ; 
  lexaddkey("Contig", &_Contig, 0) ; 
  lexaddkey("Contig9", &_Contig9, 0) ; 
  lexaddkey("Contiguous", &_Contiguous, 0) ; 
  lexaddkey("Correct_name", &_Correct_name, 0) ; 
  lexaddkey("Corresponding_DNA", &_Corresponding_DNA, 0) ; 
  lexaddkey("Corresponding_peptide", &_Corresponding_peptide, 0) ; 
  lexaddkey("Cosmid_grid", &_Cosmid_grid, 0) ; 
  lexaddkey("Cosmid_vector", &_Cosmid_vector, 0) ; 
  lexaddkey("Covers", &_Covers, 0) ; 
  lexaddkey("DB_remark", &_DB_remark, 0) ; 
  lexaddkey("DB_searched", &_DB_searched, 0) ; 
  lexaddkey("DNA", &_DNA, 0) ; 
  lexaddkey("DNA_homol", &_DNA_homol, 0) ; 
  lexaddkey("Deletion", &_Deletion, 0) ; 
  lexaddkey("Description", &_Description, 0) ; 
  lexaddkey("Designating_laboratory", &_Designating_laboratory, 0) ; 
  lexaddkey("Detection_method", &_Detection_method, 0) ; 
  lexaddkey("Df_Dup", &_Df_Dup, 0) ; 
  lexaddkey("Direct", &_Direct, 0) ; 
  lexaddkey("Display", &_Display, 0) ; 
  lexaddkey("Distance", &_Distance, 0) ; 
  lexaddkey("Does_not_hybridize_to", &_Does_not_hybridize_to, 0) ; 
  lexaddkey("Dom_let", &_Dom_let, 0) ; 
  lexaddkey("Dom_one", &_Dom_one, 0) ; 
  lexaddkey("Dom_selected", &_Dom_selected, 0) ; 
  lexaddkey("Dom_semi", &_Dom_semi, 0) ; 
  lexaddkey("Dominant", &_Dominant, 0) ; 
  lexaddkey("Duplication", &_Duplication, 0) ; 
  lexaddkey("E_mail", &_E_mail, 0) ; 
  lexaddkey("End_not_found", &_End_not_found, 0) ; 
  lexaddkey("Error", &_Error, 0) ; 
  lexaddkey("Error_scale", &_Error_scale, 0) ; 
  lexaddkey("Exact_Match_to", &_Exact_Match_to, 0) ; 
  lexaddkey("Experiment", &_Experiment, 0) ; 
  lexaddkey("Expression_construct", &_Expression_construct, 0) ; 
  lexaddkey("Extent", &_Extent, 0) ; 
  lexaddkey("Fax", &_Fax, 0) ; 
  lexaddkey("Features", &_Features, 0) ; 
  lexaddkey("FingerPrint", &_FingerPrint, 0) ; 
  lexaddkey("Finished", &_Finished, 0) ; 
  lexaddkey("Flag", &_Flag, 0) ; 
  lexaddkey("Float_Params", &_Float_Params, 0) ; 
  lexaddkey("Free_dup", &_Free_dup, 0) ; 
  lexaddkey("Freezer", &_Freezer, 0) ; 
  lexaddkey("From", &_From, 0) ; 
  lexaddkey("From_Author", &_From_Author, 0) ; 
  lexaddkey("From_Laboratory", &_From_Laboratory, 0) ; 
  lexaddkey("From_Library", &_From_Library, 0) ; 
  lexaddkey("From_left_end", &_From_left_end, 0) ; 
  lexaddkey("From_map", &_From_map, 0) ; 
  lexaddkey("Full", &_Full, 0) ; 
  lexaddkey("Full_name", &_Full_name, 0) ; 
  lexaddkey("Function", &_Function, 0) ; 
  lexaddkey("Funny_Match_to", &_Funny_Match_to, 0) ; 
  lexaddkey("Gel_Number", &_Gel_Number, 0) ; 
  lexaddkey("Gel_length", &_Gel_length, 0) ; 
  lexaddkey("Gene", &_Gene, 0) ; 
  lexaddkey("GeneId", &_GeneId, 0) ; 
  lexaddkey("Gene_class", &_Gene_class, 0) ; 
  lexaddkey("Gene_classes", &_Gene_classes, 0) ; 
  lexaddkey("General", &_General, 0) ; 
  lexaddkey("General_remark", &_General_remark, 0) ; 
  lexaddkey("Genetics", &_Genetics, 0) ; 
  lexaddkey("Genomic_Canonical", &_Genomic_Canonical, 0) ; 
  lexaddkey("Genotype", &_Genotype, 0) ; 
  lexaddkey("Gridded", &_Gridded, 0) ; 
  lexaddkey("Heat_sensitive", &_Heat_sensitive, 0) ; 
  lexaddkey("Homol", &_Homol, 0) ; 
  lexaddkey("Hybrid_cell_line", &_Hybrid_cell_line, 0) ; 
  lexaddkey("Hybridizes_to", &_Hybridizes_to, 0) ; 
  lexaddkey("In_Situ", &_In_Situ, 0) ; 
  lexaddkey("In_pool", &_In_pool, 0) ; 
  lexaddkey("Inherits", &_Inherits, 0) ; 
  lexaddkey("Inside", &_Inside, 0) ; 
  lexaddkey("Int_Params", &_Int_Params, 0) ; 
  lexaddkey("Internet", &_Internet, 0) ; 
  lexaddkey("Interval", &_Interval, 0) ; 
  lexaddkey("Intragenic_revertant_of_dominant", &_Intragenic_revertant_of_dominant, 0) ; 
  lexaddkey("Inverted", &_Inverted, 0) ; 
  lexaddkey("Is_buried_under", &_Is_buried_under, 0) ; 
  lexaddkey("Isolation", &_Isolation, 0) ; 
  lexaddkey("Isoschizomers", &_Isoschizomers, 0) ; 
  lexaddkey("Journal", &_Journal, 0) ; 
  lexaddkey("Keyword", &_Keyword, 0) ; 
  lexaddkey("Laboratory", &_Laboratory, 0) ; 
#if !defined(MACINTOSH)
  lexaddkey("Layout", &_Layout, 0) ; 
#endif
  lexaddkey("Left", &_Left, 0) ; 
  lexaddkey("Length", &_Length, 0) ; 
  lexaddkey("Library", &_Library, 0) ; 
  lexaddkey("Lines_at", &_Lines_at, 0) ; 
  lexaddkey("Linkage", &_Linkage, 0) ; 
  lexaddkey("LiquidN2", &_LiquidN2, 0) ; 
  lexaddkey("Location", &_Location, 0) ; 
  lexaddkey("Locus", &_Locus, 0) ; 
  lexaddkey("LocusId", &_LocusId, 0) ; 
  lexaddkey("Locus_1", &_Locus_1, 0) ; 
  lexaddkey("Locus_2", &_Locus_2, 0) ; 
  lexaddkey("Locus_A", &_Locus_A, 0) ; 
  lexaddkey("Locus_B", &_Locus_B, 0) ; 
  lexaddkey("Loop", &_Loop, 0) ; 
  lexaddkey("Mail", &_Mail, 0) ; 
  lexaddkey("Males", &_Males, 0) ; 
  lexaddkey("Map", &_Map, 0) ; 
  lexaddkey("Mapper", &_Mapper, 0) ; 
  lexaddkey("Mapping_data", &_Mapping_data, 0) ; 
  lexaddkey("Maps_with", &_Maps_with, 0) ; 
  lexaddkey("Marker_locus", &_Marker_locus, 0) ; 
  lexaddkey("Matching_Genomic", &_Matching_Genomic, 0) ; 
  lexaddkey("Matching_cDNA", &_Matching_cDNA, 0) ; 
  lexaddkey("Maternal", &_Maternal, 0) ; 
  lexaddkey("Max", &_Max, 0) ; 
  lexaddkey("Min", &_Min, 0) ; 
  lexaddkey("Min_score", &_Min_score, 0) ; 
  lexaddkey("Minus70", &_Minus70, 0) ; 
  lexaddkey("Molecular_information", &_Molecular_information, 0) ; 
  lexaddkey("Multi_point", &_Multi_point, 0) ; 
  lexaddkey("Multiplet", &_Multiplet, 0) ; 
  lexaddkey("Muscle", &_Muscle, 0) ; 
  lexaddkey("Mutagen", &_Mutagen, 0) ; 
  lexaddkey("N_gaps", &_N_gaps, 0) ; 
  lexaddkey("Name", &_Name, 0) ; 
  lexaddkey("Negative_clone", &_Negative_clone, 0) ; 
  lexaddkey("Negative_locus", &_Negative_locus, 0) ; 
  lexaddkey("Negative_probe", &_Negative_probe, 0) ; 
  lexaddkey("Nick_name", &_Nick_name, 0) ; 
  lexaddkey("No_overlap", &_No_overlap, 0) ; 
  lexaddkey("No_stagger", &_No_stagger, 0) ; 
  lexaddkey("Offset", &_Offset, 0) ; 
  lexaddkey("Old_CGC_distance", &_Old_CGC_distance, 0) ; 
  lexaddkey("Old_CGC_results", &_Old_CGC_results, 0) ; 
  lexaddkey("Oligo", &_Oligo, 0) ; 
  lexaddkey("One_all", &_One_all, 0) ; 
  lexaddkey("One_let", &_One_let, 0) ; 
  lexaddkey("One_recombinant", &_One_recombinant, 0) ; 
  lexaddkey("Origin", &_Origin, 0) ; 
  lexaddkey("Other_name", &_Other_name, 0) ; 
  lexaddkey("Outside", &_Outside, 0) ; 
  lexaddkey("Overhang", &_Overhang, 0) ; 
  lexaddkey("Overlap", &_Overlap, 0) ; 
  lexaddkey("PCR_product", &_PCR_product, 0) ; 
  lexaddkey("PCR_remark", &_PCR_remark, 0) ; 
  lexaddkey("Page", &_Page, 0) ; 
  lexaddkey("Paper", &_Paper, 0) ; 
  lexaddkey("Parameters", &_Parameters, 0) ; 
  lexaddkey("Paternal", &_Paternal, 0) ; 
  lexaddkey("Pattern", &_Pattern, 0) ; 
  lexaddkey("Pep_homol", &_Pep_homol, 0) ; 
  lexaddkey("Peptide", &_Peptide, 0) ; 
  lexaddkey("Percent_Identity", &_Percent_Identity, 0) ; 
  lexaddkey("Phenotype", &_Phenotype, 0) ; 
  lexaddkey("Phone", &_Phone, 0) ; 
  lexaddkey("Polymorphism", &_Polymorphism, 0) ; 
  lexaddkey("Position", &_Position, 0) ; 
  lexaddkey("Positive_clone", &_Positive_clone, 0) ; 
  lexaddkey("Positive_locus", &_Positive_locus, 0) ; 
  lexaddkey("Positive_pool_probe", &_Positive_pool_probe, 0) ; 
  lexaddkey("Positive_probe", &_Positive_probe, 0) ; 
  lexaddkey("Possible_exon", &_Possible_exon, 0) ; 
  lexaddkey("Precursor", &_Precursor, 0) ; 
  lexaddkey("Primers", &_Primers, 0) ; 
  lexaddkey("Procedure", &_Procedure, 0) ; 
  lexaddkey("Processed_mRNA", &_Processed_mRNA, 0) ; 
  lexaddkey("Product", &_Product, 0) ; 
  lexaddkey("Properties", &_Properties, 0) ; 
  lexaddkey("Protocol", &_Protocol, 0) ; 
  lexaddkey("Pseudogene", &_Pseudogene, 0) ; 
  lexaddkey("Publisher", &_Publisher, 0) ; 
  lexaddkey("Qualifier", &_Qualifier, 0) ; 
  lexaddkey("RFLP", &_RFLP, 0) ; 
  lexaddkey("RNA", &_RNA, 0) ; 
  lexaddkey("Rearrangement", &_Rearrangement, 0) ; 
  lexaddkey("Rearrangement_1", &_Rearrangement_1, 0) ; 
  lexaddkey("Rearrangement_2", &_Rearrangement_2, 0) ; 
  lexaddkey("Recessive", &_Recessive, 0) ; 
  lexaddkey("Recs_all", &_Recs_all, 0) ; 
  lexaddkey("Reference", &_Reference, 0) ; 
  lexaddkey("Reference_Allele", &_Reference_Allele, 0) ; 
  lexaddkey("Reference_strain", &_Reference_strain, 0) ; 
  lexaddkey("Refers_to", &_Refers_to, 0) ;
  lexaddkey("Rejected", &_Rejected, 0) ;
  lexaddkey("Rejected_from", &_Rejected_from, 0) ;
  lexaddkey("Rejected_Reads", &_Rejected_Reads, 0) ;
  lexaddkey("Related_Sequence", &_Related_Sequence, 0) ; 
  lexaddkey("Remark", &_Remark, 0) ; 
  lexaddkey("Repeats", &_Repeats, 0) ; 
  lexaddkey("Representative", &_Representative, 0) ; 
  lexaddkey("Results", &_Results, 0) ; 
  lexaddkey("Right", &_Right, 0) ; 
  lexaddkey("Row", &_Row, 0) ; 
  lexaddkey("Score", &_Score, 0) ; 
  lexaddkey("Search_Method", &_Search_Method, 0) ; 
  lexaddkey("Selected", &_Selected, 0) ; 
  lexaddkey("Selected_loci", &_Selected_loci, 0) ; 
  lexaddkey("Selected_trans", &_Selected_trans, 0) ; 
  lexaddkey("Semi_dominant", &_Semi_dominant, 0) ; 
  lexaddkey("Seq_length", &_Seq_length, 0) ; 
  lexaddkey("Sequence", &_Sequence, 0) ; 
  lexaddkey("Sequencing_vector", &_Sequencing_vector, 0) ; 
  lexaddkey("Sex_cis", &_Sex_cis, 0) ; 
  lexaddkey("Sex_full", &_Sex_full, 0) ; 
  lexaddkey("Sex_one", &_Sex_one, 0) ; 
  lexaddkey("Shotgun", &_Shotgun, 0) ; 
  lexaddkey("Significant_bases", &_Significant_bases, 0) ; 
  lexaddkey("Similarity", &_Similarity, 0) ; 
  lexaddkey("Site", &_Site, 0) ; 
  lexaddkey("Source", &_Source, 0) ; 
  lexaddkey("Source_Exons", &_Source_Exons, 0) ; 
  lexaddkey("Space_at", &_Space_at, 0) ; 
  lexaddkey("Species", &_Species, 0) ; 
  lexaddkey("Staff", &_Staff, 0) ; 
  lexaddkey("Start_not_found", &_Start_not_found, 0) ; 
#if !defined(MACINTOSH)
  lexaddkey("Status", &_Status, 0) ; 
#endif
  lexaddkey("Stop", &_Stop, 0) ; 
  lexaddkey("Strain", &_Strain, 0) ; 
  lexaddkey("Strain_designation", &_Strain_designation, 0) ; 
  lexaddkey("Strictly_Maternal", &_Strictly_Maternal, 0) ; 
  lexaddkey("Structure", &_Structure, 0) ; 
  lexaddkey("Subsequence", &_Subsequence, 0) ; 
  lexaddkey("TATA_signal", &_TATA_signal, 0) ; 
  lexaddkey("TSL", &_TSL, 0) ; 
  lexaddkey("TSL_site", &_TSL_site, 0) ; 
  lexaddkey("Tandem", &_Tandem, 0) ; 
  lexaddkey("Temperature", &_Temperature, 0) ; 
  lexaddkey("Temperature_sensitive", &_Temperature_sensitive, 0) ; 
  lexaddkey("Tested", &_Tested, 0) ; 
  lexaddkey("Text_Params", &_Text_Params, 0) ; 
  lexaddkey("Title", &_Title, 0) ; 
  lexaddkey("To_right_end", &_To_right_end, 0) ; 
  lexaddkey("Trace_quality", &_Trace_quality, 0) ; 
  lexaddkey("Excellent_upto", &_Excellent_upto, 0) ; 
  lexaddkey("Good_upto", &_Good_upto, 0) ; 
  lexaddkey("Fair_upto", &_Fair_upto, 0) ; 
  lexaddkey("Transcript", &_Transcript, 0) ; 
  lexaddkey("Translocation", &_Translocation, 0) ; 
  lexaddkey("Transposon", &_Transposon, 0) ; 
  lexaddkey("Transposon_insertion", &_Transposon_insertion, 0) ; 
  lexaddkey("Type", &_Type, 0) ; 
  lexaddkey("Type_1", &_Type_1, 0) ; 
  lexaddkey("Type_2", &_Type_2, 0) ; 
  lexaddkey("Unit_Length", &_Unit_Length, 0) ; 
  lexaddkey("Unprocessed_mRNA", &_Unprocessed_mRNA, 0) ; 
  lexaddkey("Variant", &_Variant, 0) ; 
  lexaddkey("Variant_of", &_Variant_of, 0) ; 
  lexaddkey("Vaxmap", &_Vaxmap, 0) ; 
  lexaddkey("Vector", &_Vector, 0) ; 
  lexaddkey("Version", &_Version, 0) ; 
  lexaddkey("Volume", &_Volume, 0) ; 
  lexaddkey("Weak", &_Weak, 0) ; 
  lexaddkey("Well_ordered", &_Well_ordered, 0) ; 
  lexaddkey("With_Maternal_Effect", &_With_Maternal_Effect, 0) ; 
  lexaddkey("Y_remark", &_Y_remark, 0) ; 
  lexaddkey("Year", &_Year, 0) ; 
  lexaddkey("cDNA", &_cDNA, 0) ; 
  lexaddkey("mRNA", &_mRNA, 0) ; 
  lexaddkey("Mat_peptide", &_mat_peptide, 0) ; 
  lexaddkey("Misc_feature", &_misc_feature, 0) ; 
  lexaddkey("Misc_signal", &_misc_signal, 0) ; 
  lexaddkey("Modified_base", &_modified_base, 0) ; 
  lexaddkey("Mutation", &_mutation, 0) ; 
  lexaddkey("Old_sequence", &_old_sequence, 0) ; 
  lexaddkey("pMap", &_pMap, 0) ; 
  lexaddkey("polyA_signal", &_polyA_signal, 0) ; 
  lexaddkey("polyA_site", &_polyA_site, 0) ; 
  lexaddkey("Promoter", &_promoter, 0) ; 
  lexaddkey("rRNA", &_rRNA, 0) ; 
  lexaddkey("Repeat_region", &_repeat_region, 0) ; 
  lexaddkey("Repeat_unit", &_repeat_unit, 0) ; 
  lexaddkey("Sig_peptide", &_sig_peptide, 0) ; 
  lexaddkey("snRNA", &_snRNA, 0) ; 
  lexaddkey("tRNA", &_tRNA, 0) ; 
  lexaddkey("A_Repeat", &_A_Repeat, 0) ; 
  lexaddkey("Note", &_Note, 0) ; 
  lexaddkey("Subpool", &_Subpool, 0) ; 
  lexaddkey("Colour", &_Colour, 0) ; 
  lexaddkey("ORF", &_ORF, 0) ; 
  lexaddkey("ATG", &_ATG, 0) ; 
  lexaddkey("Splice3", &_Splice3, 0) ; 
  lexaddkey("Splice5", &_Splice5, 0) ; 
  lexaddkey("Coding_seg", &_Coding_seg, 0) ; 
  lexaddkey("Foreign_Reference", &_Foreign_Reference, 0) ; 
  lexaddkey("Virtual_row", &_Virtual_row, 0) ; 
  lexaddkey("In_grid", &_In_grid, 0) ; 
  lexaddkey("Chrom_Band", &_Chrom_Band, 0) ; 
  lexaddkey("Centromere", &_Centromere, 0) ; 
  lexaddkey("Dark", &_Dark, 0) ; 
  lexaddkey("NOR", &_NOR, 0) ; 
  lexaddkey("Inherits_from", &_Inherits_from, 0) ; 
  lexaddkey("Fate", &_Fate, 0) ; 
  lexaddkey("Lineage", &_Lineage, 0) ; 
  lexaddkey("Parent", &_Parent, 0) ; 
  lexaddkey("Daughter", &_Daughter, 0) ; 
  lexaddkey("Group", &_Group, 0) ; 
  lexaddkey("In_group", &_In_group, 0) ; 
  lexaddkey("Group_member", &_Group_member, 0) ; 
  lexaddkey("Neurodata", &_Neurodata, 0) ; 
  lexaddkey("Send", &_Send, 0) ; 
  lexaddkey("Send_joint", &_Send_joint, 0) ; 
  lexaddkey("Receive", &_Receive, 0) ; 
  lexaddkey("Receive_joint", &_Receive_joint, 0) ; 
  lexaddkey("Gap_junction", &_Gap_junction, 0) ; 
  lexaddkey("Contact", &_Contact, 0) ; 
  lexaddkey("Lineage_name", &_Lineage_name, 0) ; 
  lexaddkey("Equivalence_fate", &_Equivalence_fate, 0) ; 
  lexaddkey("Equivalence_origin", &_Equivalence_origin, 0) ; 
  lexaddkey("Embryo_division_time", &_Embryo_division_time, 0) ; 
  lexaddkey("Clone_left_end", &_Clone_left_end, 0) ; 
  lexaddkey("Clone_right_end", &_Clone_right_end, 0) ; 
  lexaddkey("Derived_from", &_Derived_from, 0) ; 
  lexaddkey("Derivative", &_Derivative, 0) ; 
  lexaddkey("Isolated_for", &_Isolated_for, 0) ; 
  lexaddkey("Wormpep", &_Wormpep, 0) ; 
  lexaddkey("Replaces", &_Replaces, 0) ; 
  lexaddkey("Replaced_by", &_Replaced_by, 0) ; 
  lexaddkey("Confidential_remark", &_Confidential_remark, 0) ; 
  lexaddkey("Submitted", &_Submitted, 0) ; 
  lexaddkey("Overlap_right", &_Overlap_right, 0) ; 
  lexaddkey("Overlap_left", &_Overlap_left, 0) ; 
  lexaddkey("Repeat_consensus", &_Repeat_consensus, 0) ; 
  lexaddkey("cDNA_EST", &_cDNA_EST, 0) ; 
  lexaddkey("Archived", &_Archived, 0) ; 
  lexaddkey("Lethal", &_Lethal, 0) ; 
  lexaddkey("Hand_fixed", &_Hand_fixed, 0) ; 
  lexaddkey("Balancer", &_Balancer, 0) ; 
  lexaddkey("Region", &_Region, 0) ; 
  lexaddkey("Hand_verified", &_Hand_verified, 0) ; 
  lexaddkey("YAC", &_YAC, 0) ; 
  lexaddkey("Cytogenetic", &_Cytogenetic, 0) ; 
  lexaddkey("Inside_YAC", &_Inside_YAC, 0) ; 
  lexaddkey("Inside_Fragment", &_Inside_Fragment, 0) ; 
  lexaddkey("Bibliography", &_Bibliography, 0) ; 
  lexaddkey("Links", &_Links, 0) ; 
  lexaddkey("Link", &_Link, 0) ; 
  lexaddkey("Chimeric", &_Chimeric, 0) ; 
  lexaddkey("Non_Chimeric", &_Non_Chimeric, 0) ; 
  lexaddkey("Size", &_Size, 0) ; 
  lexaddkey("Interval_Mapping", &_Interval_Mapping, 0) ; 
  lexaddkey("Main_Marker", &_Main_Marker, 0) ; 
  lexaddkey("Ends", &_Ends, 0) ; 
  lexaddkey("Mapping", &_Mapping, 0) ; 
  lexaddkey("p_Telomere", &_p_Telomere, 0) ; 
  lexaddkey("q_Telomere", &_q_Telomere, 0) ; 
  lexaddkey("Drawing", &_Drawing, 0) ; 
  lexaddkey("Multi_Position", &_Multi_Position, 0) ; 
  lexaddkey("Does_not_Contain", &_Does_not_Contain, 0) ; 
  lexaddkey("Contained_Locus", &_Contained_Locus, 0) ; 
  lexaddkey("Exterior_Locus", &_Exterior_Locus, 0) ; 
  lexaddkey("Inside_Rearr", &_Inside_Rearr, 0) ; 
  lexaddkey("Outside_Rearr", &_Outside_Rearr, 0) ; 
  lexaddkey("Overlaps", &_Overlaps, 0) ; 
  lexaddkey("Does_not_Overlap", &_Does_not_Overlap, 0) ; 
  lexaddkey("No_Overlap_Rearr", &_No_Overlap_Rearr, 0) ; 
  lexaddkey("Overlapping_Rearr", &_Overlapping_Rearr, 0) ; 
  lexaddkey("Contained_Rearr", &_Contained_Rearr, 0) ; 
  lexaddkey("Email", &_Email, 0) ; 
  lexaddkey("Flipped", &_Flipped, 0) ;
  lexaddkey("Assembled_as_dummy", &_Previous_contig, 0) ;
  lexaddkey("Previous_contig", &_Previous_contig, 0) ;
  lexaddkey("Assembled_from", &_Assembled_from, 0) ; 
  lexaddkey("Assembled_into", &_Assembled_into, 0) ; 
  lexaddkey("Proposed", &_Proposed, 0) ; 
  lexaddkey("BaseCall", &_BaseCall, 0) ; 
  lexaddkey("Genetic", &_Genetic, 0) ; 
  lexaddkey("Physical", &_Physical, 0) ; 
  lexaddkey("Restriction_Enzyme", &_Restriction_Enzyme, 0) ; 
  lexaddkey("AA", &_AA, 0) ; 
  lexaddkey("Clone_Grid", &_Clone_Grid, 0) ; 
  lexaddkey("Df_Dup_data", &_Df_Dup_data, 0) ; 
  lexaddkey("Fragment", &_Fragment, 0) ; 
  lexaddkey("Pages", &_Pages, 0) ; 
  lexaddkey("Medline_ID", &_Medline_ID, 0) ; 
  lexaddkey("Lab_Location", &_Lab_Location, 0) ; 
  lexaddkey("Match", &_Match, 0) ; 
  lexaddkey("Method", &_Method, 0) ; 
  lexaddkey("Motif", &_Motif, 0) ; 
  lexaddkey("MultiMap", &_MultiMap, 0) ; 
  lexaddkey("Multi_counts", &_Multi_counts, 0) ; 
  lexaddkey("Multi_pt_data", &_Multi_pt_data, 0) ; 
  lexaddkey("Pool", &_Pool, 0) ; 
  lexaddkey("Unpublished", &_Unpublished, 0) ; 
  lexaddkey("EC_enzyme", &_EC_enzyme, 0) ; 
  lexaddkey("Document", &_Document, 0) ; 
  lexaddkey("Unit", &_Unit, 0) ; 
  lexaddkey("Enzyme", &_Enzyme, 0) ; 
  lexaddkey("GDB_id", &_GDB_id, 0) ; 
  lexaddkey("OMIM", &_OMIM, 0) ; 
  lexaddkey("Probe", &_Probe, 0) ; 
  lexaddkey("ATCC_ref", &_ATCC_ref, 0) ; 
  lexaddkey("Add_date", &_Add_date, 0) ; 
  lexaddkey("Annotation", &_Annotation, 0) ; 
  lexaddkey("Approval_date", &_Approval_date, 0) ; 
  lexaddkey("Approved", &_Approved, 0) ; 
  lexaddkey("Assign_mode", &_Assign_mode, 0) ; 
  lexaddkey("Availability", &_Availability, 0) ; 
  lexaddkey("CS", &_CS, 0) ; 
  lexaddkey("Cosmid_clones", &_Cosmid_clones, 0) ; 
  lexaddkey("Created", &_Created, 0) ; 
  lexaddkey("Current_symbol", &_Current_symbol, 0) ; 
  lexaddkey("Cytogenetic_location", &_Cytogenetic_location, 0) ; 
  lexaddkey("DNA_type", &_DNA_type, 0) ; 
  lexaddkey("Data_source", &_Data_source, 0) ; 
  lexaddkey("EC_ref", &_EC_ref, 0) ; 
  lexaddkey("Editing_Date", &_Editing_Date, 0) ; 
  lexaddkey("Excision_Data", &_Excision_Data, 0) ; 
  lexaddkey("GDB", &_GDB, 0) ; 
  lexaddkey("GenBank_ref", &_GenBank_ref, 0) ; 
  lexaddkey("HGML_ref", &_HGML_ref, 0) ; 
  lexaddkey("KW", &_KW, 0) ; 
  lexaddkey("Links_to", &_Links_to, 0) ; 
  lexaddkey("Max_heterozygosity", &_Max_heterozygosity, 0) ; 
  lexaddkey("Modification_date", &_Modification_date, 0) ; 
  lexaddkey("Modified", &_Modified, 0) ; 
  lexaddkey("OMIM_document", &_OMIM_document, 0) ; 
  lexaddkey("Old_Num", &_Old_Num, 0) ; 
  lexaddkey("Originator", &_Originator, 0) ; 
  lexaddkey("Polymorphism_Data", &_Polymorphism_Data, 0) ; 
  lexaddkey("Previous_symbol", &_Previous_symbol, 0) ; 
  lexaddkey("Vector_type", &_Vector_type, 0) ; 
  lexaddkey("Primer", &_Primer, 0) ; 
  lexaddkey("SCF_File", &_SCF_File, 0) ; 
  lexaddkey("PCR", &_PCR, 0) ; 
  lexaddkey("Cor", &_Cor, 0) ; 
  lexaddkey("Neighbours", &_Neighbours, 0) ; 
  lexaddkey("Common_bands", &_Common_bands, 0) ; 
  lexaddkey("Band_Values", &_Band_Values, 0) ; 
  lexaddkey("Segment_Lengths", &_Segment_Lengths, 0) ; 
  lexaddkey("Fragment_of", &_Fragment_of, 0) ; 
  lexaddkey("ABI", &_ABI, 0) ; 
  lexaddkey("Wild_type", &_Wild_type, 0) ; 
  lexaddkey("ABI_Date", &_ABI_Date, 0) ; 
  lexaddkey("ABI_Comment", &_ABI_Comment, 0) ; 
  lexaddkey("ABI_Machine", &_ABI_Machine, 0) ; 
  lexaddkey("Sample", &_Sample, 0) ; 
  lexaddkey("ABI_Analysis", &_ABI_Analysis, 0) ; 
  lexaddkey("Run_Start", &_Run_Start, 0) ; 
  lexaddkey("Run_Stop", &_Run_Stop, 0) ; 
  lexaddkey("Clipped", &_Clipped, 0) ; 
  lexaddkey("Old_Clipping", &_Old_Clipping, 0) ; 
  lexaddkey("Hand_Clipping", &_Hand_Clipping, 0) ; 
  lexaddkey("Vector_Clipping", &_Vector_Clipping, 0) ; 
  lexaddkey("Clipping", &_Clipping, 0) ; 
  lexaddkey("Symbol", &_Symbol, 0) ; 
  lexaddkey("Anchor", &_Anchor, 0) ; 
  lexaddkey("Band", &_Band, 0) ; 
  lexaddkey("Gel", &_Gel, 0) ; 
  lexaddkey("Band_Lengths", &_Band_Lengths, 0) ; 
  lexaddkey("Otto", &_Otto, 0) ; 
  lexaddkey("Main_Locus", &_Main_Locus, 0) ; 
  lexaddkey("STS", &_STS, 0) ; 
  lexaddkey("SELF", &_SELF, 0) ; 
  lexaddkey("Checked_individually", &_Checked_individually, 0) ; 
  lexaddkey("Outside_YAC", &_Outside_YAC, 0) ; 
  lexaddkey("STS_out", &_STS_out, 0) ; 
  lexaddkey("AluPCR", &_AluPCR, 0) ; 
  lexaddkey("by_AluPCR_Hybridisation", &_by_AluPCR_Hybridisation, 0) ; 
  lexaddkey("by_Finger_Printing", &_by_Finger_Printing, 0) ; 
  lexaddkey("Kpn", &_Kpn, 0) ; 
  lexaddkey("THE", &_THE, 0) ; 
  lexaddkey("Somatic_hybrid", &_Somatic_hybrid, 0) ; 
  lexaddkey("Genethon", &_Genethon, 0) ; 
  lexaddkey("Route_2", &_Route_2, 0) ; 
  lexaddkey("Inside_YAC_Not_Verified", &_Inside_YAC_Not_Verified, 0) ; 
  lexaddkey("Locus_Not_Verified", &_Locus_Not_Verified, 0) ; 
  lexaddkey("Possibly_Contains_Locus", &_Possibly_Contains_Locus, 0) ; 
  lexaddkey("Possibly_Inside_YAC", &_Possibly_Inside_YAC, 0) ; 
  lexaddkey("Locus_Out", &_Locus_Out, 0) ; 
  lexaddkey("Allelic_Variants", &_Allelic_Variants, 0) ; 
  lexaddkey("Clinical_Synopses", &_Clinical_Synopses, 0) ; 
  lexaddkey("Cross_References", &_Cross_References, 0) ; 
  lexaddkey("Diagnosis", &_Diagnosis, 0) ; 
  lexaddkey("Evolution", &_Evolution, 0) ; 
  lexaddkey("Historical_Information", &_Historical_Information, 0) ; 
  lexaddkey("Index_Terms", &_Index_Terms, 0) ; 
  lexaddkey("Last_Edited", &_Last_Edited, 0) ; 
  lexaddkey("Locus_Description", &_Locus_Description, 0) ; 
  lexaddkey("MIM_Number", &_MIM_Number, 0) ; 
  lexaddkey("Old_MIM_Number", &_Old_MIM_Number, 0) ; 
  lexaddkey("Population_Genetics", &_Population_Genetics, 0) ; 
  lexaddkey("References", &_References, 0) ; 
  lexaddkey("Treatment", &_Treatment, 0) ; 
  lexaddkey("Main_Text", &_Main_Text, 0) ; 
  lexaddkey("Animal_Models", &_Animal_Models, 0) ; 
  lexaddkey("Map_location", &_map_location, 0) ; 
  lexaddkey("Map_error", &_map_error, 0) ; 
  lexaddkey("LongText", &_LongText, 0) ; 
  lexaddkey("PCR_data", &_PCR_data, 0) ; 
  lexaddkey("Clone_data", &_Clone_data, 0) ; 
  lexaddkey("Oligo_DNA", &_Oligo_DNA, 0) ; 
  lexaddkey("Allelep", &_Allelep, 0) ; 
  lexaddkey("Allelem", &_Allelem, 0) ; 
  lexaddkey("Gametep", &_Gametep, 0) ; 
  lexaddkey("Gametem", &_Gametem, 0) ; 
  lexaddkey("Tetrad", &_Tetrad, 0) ; 
  lexaddkey("Centromere_segregation", &_Centromere_segregation, 0) ; 
  lexaddkey("Genetic_code",  &_Genetic_code, 0) ;
  lexaddkey("Lethal_tested", &_Lethal_tested, 0) ; 

  lexaddkey("Description", &_description, 0) ;
  lexaddkey("Used_in", &_used_in, 0) ;
  lexaddkey("Action_type", &_action_type, 0) ;
  lexaddkey("Internal", &_internal, 0) ;
  lexaddkey("Ace_query", &_ace_query, 0) ;
  lexaddkey("Ace_dump", &_ace_dump, 0) ;
  lexaddkey("Name_dump", &_name_dump, 0) ;
  lexaddkey("Tace_command", &_tace_command, 0) ;
  lexaddkey("External", &_external, 0) ;
  lexaddkey("Compound", &_compound, 0) ;
  lexaddkey("Action_element", &_action_element, 0) ;
  lexaddkey("Params", &_params, 0) ;
  lexaddkey("No_report_status", &_no_report_status, 0) ;
  lexaddkey("Read_ace", &_read_ace, 0) ;
  lexaddkey("Synchro", &_synchro, 0) ;
  lexaddkey("Import_keyset", &_import_keyset, 0) ;

/* for comparative maps (added by jld) */

  lexaddkey ("Homology", &_Homology, 0) ;
  lexaddkey ("Doc", &_Doc, 0) ;
  lexaddkey ("Pairwise", &_Pairwise, 0) ;
  lexaddkey ("Fuzzy", &_Fuzzy, 0) ;
  lexaddkey ("Symbol", &_Symbol, 0) ;

  lexaddkey ("SMAP", &_SMAP, 0) ;
  lexaddkey ("S_Parent", &_S_Parent, 0) ;
  /* lexaddkey ("", &_, 0) ; */
}

/***********************************************************************************/


int
  _VPeptide,
  _VSequence,
  _VProtein,
  _VDNA,
  _VPaper,
  _VMethod,

  _VMap,
  _VgMap,
  _VvMap,
  _VMultiMap,

  _VLocus,
  _VGene,
  _VAllele,
  _VInterval,
  _V2_point_data,
  _VMulti_pt_data,

  _VClone,
  _VClone_Grid,
  _VPool,
  _VContig,
  _VpMap,

  _VChrom_Band,
  _VMotif,
  _VBaseCall,
  _VBaseQuality,
  _VBasePosition,
  _VOligoRepeat,
  _VPfam,

/* for comparative maps (added by jld) */

  _VHomology_group,
  _VMap_set,
  _VDoc,
  _VGenetic_code ,
  _VPerson ;

void classInit (void)
{
  KEY key ;
  extern BOOL READING_MODELS ;
  BOOL old = READING_MODELS ;
  READING_MODELS = TRUE ;
 
  lexaddkey ("Peptide", &key, _VMainClasses) ; _VPeptide = KEYKEY(key) ;
  lexaddkey ("Sequence", &key, _VMainClasses) ; _VSequence = KEYKEY(key) ;
  lexaddkey ("Protein", &key, _VMainClasses) ; _VProtein = KEYKEY(key) ;
  lexaddkey ("DNA", &key, _VMainClasses) ; _VDNA = KEYKEY(key) ;
  lexaddkey ("Paper", &key, _VMainClasses) ; _VPaper = KEYKEY(key) ;
  lexaddkey ("Method", &key, _VMainClasses) ; _VMethod = KEYKEY(key) ;
 
  lexaddkey ("Map", &key, _VMainClasses) ; _VMap = KEYKEY(key) ;
  lexaddkey ("gMap", &key, _VMainClasses) ; _VgMap = KEYKEY(key) ;
  lexaddkey ("vMap", &key, _VMainClasses) ; _VvMap = KEYKEY(key) ;
  lexaddkey ("MultiMap", &key, _VMainClasses) ; _VMultiMap = KEYKEY(key) ;

  lexaddkey ("Locus", &key, _VMainClasses) ; _VLocus = KEYKEY(key) ;
  lexaddkey ("Gene", &key, _VMainClasses) ; _VGene = KEYKEY(key) ;
  lexaddkey ("Allele", &key, _VMainClasses) ; _VAllele = KEYKEY(key) ;
  lexaddkey ("Interval", &key, _VMainClasses) ; _VInterval = KEYKEY(key) ;
  lexaddkey ("2_point_data", &key, _VMainClasses) ; _V2_point_data = KEYKEY(key) ;
  lexaddkey ("Multi_pt_data", &key, _VMainClasses) ; _VMulti_pt_data = KEYKEY(key) ;

  lexaddkey ("Clone", &key, _VMainClasses) ; _VClone = KEYKEY(key) ;
  lexaddkey ("Clone_Grid", &key, _VMainClasses) ; _VClone_Grid = KEYKEY(key) ;
  lexaddkey ("Pool", &key, _VMainClasses) ; _VPool = KEYKEY(key) ;
  lexaddkey ("Contig", &key, _VMainClasses) ; _VContig = KEYKEY(key) ;
  lexaddkey ("pMap", &key, _VMainClasses) ; _VpMap = KEYKEY(key) ;

  lexaddkey ("Chrom_Band", &key, _VMainClasses) ; _VChrom_Band = KEYKEY(key) ;
  lexaddkey ("Motif", &key, _VMainClasses) ; _VMotif = KEYKEY(key) ;
  lexaddkey ("BaseCall", &key, _VMainClasses) ; _VBaseCall = KEYKEY(key) ;
  lexaddkey ("BaseQuality", &key, _VMainClasses) ; _VBaseQuality = KEYKEY(key) ;
  lexaddkey ("BasePosition", &key, _VMainClasses) ; _VBasePosition = KEYKEY(key) ;
  lexaddkey ("OligoRepeat", &key, _VMainClasses) ; _VOligoRepeat = KEYKEY(key) ;
  lexaddkey ("Pfam", &key, _VMainClasses) ; _VPfam = KEYKEY(key) ;

/* for comparative maps (added by jld) */

  lexaddkey ("Homology_group", &key, _VMainClasses) ; _VHomology_group = KEYKEY (key) ;
  lexaddkey ("Map_set", &key, _VMainClasses) ; _VMap_set = KEYKEY (key) ;
  lexaddkey ("Doc", &key, _VMainClasses) ; _VDoc = KEYKEY (key) ;
  lexaddkey ("Genetic_code", &key, _VMainClasses) ; _VGenetic_code = KEYKEY (key) ;
  lexaddkey ("Person", &key, _VMainClasses) ; _VPerson = KEYKEY (key) ;

  READING_MODELS = old ;

  return ;
}

/********************* end of file *********************/
 
 
