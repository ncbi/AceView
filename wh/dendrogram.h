/*  dendrogram.h : header file for dendrogram.c
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Defines the class model behaviors for new ?Tree class
 *		which can draw dendrograms of phylogenetic or taxonomic trees.
 *
 * Exported:
 *	void readTreeFile(int treeType) // reads in a New Hampshire formatted tree
 *	void readTaxonomicTree(void) ; // Specialized versions of the above...
 *	void readDNATree(void);
 *	void readProteinTree(void) ;
 *
 * HISTORY:
 * Last edited: Aug 31, 1998 (rbrusk): CELL_LINEAGE tree type added 
 * Created: Jun 9 15:22 1998 (rbrusk)
 *-------------------------------------------------------------------
 */
  
/* Some possible treeType -- one can add to these, 
   but one should modify the code accordingly */
#define UNKNOWN_TREE    0
#define TAXONOMIC_TREE  1
#define DNA_TREE        2
#define PROTEIN_TREE    4
#define CELL_LINEAGE    8

void readTreeFile(int treeType) ; /* treeType is one of the above */
void readTaxonomicTree(void) ;   
void readDNATree(void);   
void readProteinTree(void) ;   
