
#include "../wac/ac.h"
#include "vtxt.h"
#include "vtxt_.h"
#include "utils.h"

/*************************************************************/

static void geneXmlExportXmltag2Acetag (char **cpp0, vTXT buf, char *xmlTag, char *aceTag, BOOL multi, AC_HANDLE h)
{
  char *data ;

/* usage: cp = xml ; while (cq = xmlGetNextTagContent (&cp, tag, h)) print (cq) ; */
/* extracts iterativelly the content of all occurences of tag in xml */

  while ((data =  xmlGetNextTagContent (cpp0, xmlTag, h)))
    {
      vtxtPrintf (buf, "%s %s\n", aceTag, freeprotect (data)) ; 
      if (!multi) break ;
    }
} /*  geneXmlExportXmltag2Acetag */

/*************************************************************/
/*
<Gene-commentary_type value="phenotype">19</Gene-commentary_type>
<Gene-commentary_text>LADD syndrome</Gene-commentary_text>
<Gene-commentary>
<Gene-commentary_type value="phenotype">19</Gene-commentary_type>
<Gene-commentary_text>Bladder cancer</Gene-commentary_text>
<Gene-commentary_version>0</Gene-commentary_version>
<Gene-commentary_source>
<Other-source>
<Other-source_src>
<Dbtag>
<Dbtag_db>MIM</Dbtag_db>
<Dbtag_tag>
<Object-id>
<Object-id_id>109800</Object-id_id>
</Object-id>
</Dbtag_tag>
</Dbtag>
</Other-source_src>
<Other-source_anchor>MIM: 109800</Other-source_anchor>
</Other-source>
</Gene-commentary_source>

  Locus_Phenotype     "[OMIM 109800] Bladder cancer"

 KEGG_pathway "MAPK signaling pathway"   -> function dans le web


<Gene-commentary>
<Gene-commentary_type value="comment">254</Gene-commentary_type>
<Gene-commentary_heading>Pathways</Gene-commentary_heading>
<Gene-commentary_version>0</Gene-commentary_version>
<Gene-commentary_comment>
<Gene-commentary>
<Gene-commentary_type value="comment">254</Gene-commentary_type>
<Gene-commentary_text>KEGG pathway: MAPK signaling pathway</Gene-commentary_text>
<Gene-commentary_version>0</Gene-commentary_version>
<Gene-commentary_source>
<Other-source>
<Other-source_src>
<Dbtag>
<Dbtag_db>KEGG pathway</Dbtag_db>
<Dbtag_tag>
<Object-id>
<Object-id_str>04010</Object-id_str>
</Object-id>
</Dbtag_tag>
</Dbtag>
</Other-source_src>
<Other-source_anchor>04010</Other-source_anchor>
<Other-source_url>http://www.genome.jp/dbget-bin/show_pathway?hsa04010+2261</Other-source_url>
</Other-source>
</Gene-commentary_source>
</Gene-commentary>

*/

/*************************************************************/

static void geneXmlExportKegg (const char *geneid, char **cpp0, vTXT buf, AC_HANDLE h)
{
  char *block, *cp, *data, *kegg ;
  
  while ((block =  xmlGetNextTagContent (cpp0, "Gene-commentary", h)))
    {
      cp = block ;
      while ((data = xmlGetNextTagContent (&cp, "Gene-commentary_text", h)))
	{
	  if (strstr (data, "KEGG pathway"))
	    {
	      data += strlen ("KEGG pathway") + 2 ;
	      /* cp = block ; */
	      if ((kegg = xmlGetNextTagContent (&cp, "Object-id_str", h)))
		{
		  vtxtPrintf (buf, "\nExtern KEGG_%s\nKEGG\nPathway %s\n\n"
			      , kegg, freeprotect (data)) ;
		  vtxtPrintf (buf, "GeneId %s\nExtern KEGG_%s\nPathway %s\n", geneid, kegg, freeprotect (data)) ;
		}
	    }
	}
    }
} /* geneXmlExportKegg */

/*************************************************************/

static void geneXmlExportAccession (char **cpp0, vTXT buf, AC_HANDLE h)
{
  char *block, *acc ;
  
  while ((block =  xmlGetNextTagContent (cpp0, "Gene-commentary", h)))
    {
      if (xmlGetNextTagContent (&block, "Gene-commentary_type value", h) &&
	  (acc = xmlGetNextTagContent (&block, "Gene-commentary_accession", h)) &&
	  (
	   !strncmp (acc, "NM_",3) ||
	   !strncmp (acc, "NR_",3) ||
	   !strncmp (acc, "XM_",3) ||
	   !strncmp (acc, "XR_",3)
	   )
	  )
	vtxtPrintf (buf, "Sequence %s\n", freeprotect (acc)) ;
    }
} /* geneXmlExportAccession */

/*************************************************************/

/* parsing:
 *  <Gene-commentary_label>EC</Gene-commentary_label><Gene-commentary_text>xxxx<</Gene-commentary_text>
 */
static void geneXmlExportXmlTagValue2Acetag (char **cpp0, vTXT buf
					     , char *xmlTag1, char *xmlTag2
					     , char *aceTag, AC_HANDLE h)
{
  char *block, *tag ;
  while ((block =  xmlGetNextTagContent (cpp0, xmlTag1, h)))
    {
      while ((tag =  xmlGetNextTagContent (&block, messprintf ("%s_label", xmlTag1), h)))
	if (! strcmp (tag, xmlTag2))
	  geneXmlExportXmltag2Acetag (&block, buf, messprintf ("%s_text", xmlTag1), aceTag, FALSE, h) ;
    }
} /*  geneXmlExportXmlTagValue2Acetag */

/*************************************************************/

/* parsing:
 * <Gene-commentary_type value="phenotype">19</Gene-commentary_type>
 * <Gene-commentary_text>LADD syndrome</Gene-commentary_text>
 * <Gene-commentary>
 */
static void geneXmlExportXmlTagWithValue2Acetag (const char *geneid 
						 , char **cpp0, vTXT buf
						 , char *xmlTag1, char *xmlTag2
						 , char *aceTag, AC_HANDLE h)
{
  char *block, *cp0, *origin, *omim, *data ;
  while ((block =  xmlGetNextTagContent (cpp0, xmlTag1, h)))
    {
      if (xmlGetNextTagContent (&block, messprintf ("%s_type value=\"%s\"", xmlTag1, xmlTag2), h))
	{
	  cp0 = block ; omim = 0 ;
	  origin = xmlGetNextTagContent (&cp0, "Dbtag_db", h) ;
	  if (origin)
	    omim =  xmlGetNextTagContent (&cp0, "Object-id_id", h) ;  
	  if  ((data =  xmlGetNextTagContent (&block, messprintf ("%s_text", xmlTag1), h)))
	    {
	      if (omim && strstr (origin, "MIM"))
		{
		  vtxtPrintf (buf, "\nExtern  \"OMIM_%s\"\nOMIM\nOMIM_alias %s %s\n\n"
			      , omim, freeprotect (data), geneid
			      ) ;
		  vtxtPrintf (buf, "geneId %s\n", geneid) ; /* go back to geneid */
		  vtxtPrintf (buf, "Extern \"OMIM_%s\"\n", omim) ;
		  vtxtPrintf (buf, "Disease %s\n",  freeprotect (data)) ;
		}
 	      else if (omim)  
		vtxtPrintf (buf, "%s %s\n", aceTag
			    , freeprotect (messprintf ("[%s %s] %s", origin, omim, data))) ; 
	      else 
		vtxtPrintf (buf, "%s %s\n", aceTag, freeprotect (data)) ; 
	    }
	}
    }
} /*  geneXmlExportXmlTagValue2Acetag */

/*************************************************************/

static char *geneXmlGeneContent (const char *geneid, char *cp0, vTXT buf, int line, AC_HANDLE h)
{
  char *cp, *cq, *summary, *geneName, *generef ;

  cp = cp0 ; geneName = xmlGetNextTagContent (&cp, "Gene-ref_locus", h) ;
  if (geneName)
    {
      vtxtPrintf (buf, "LocusLink \"%s\"\n", geneName) ;
      cp = cp0 ; summary = xmlGetNextTagContent (&cp, "Entrezgene_summary", h) ;
      if (summary)
	vtxtPrintf (buf, "Summary %s\n"
		    , freeprotect (messprintf ("[RefSeq Summary %s] %s", geneName ? geneName : "", summary))) ; 
      cp = cp0 ; geneXmlExportXmlTagValue2Acetag (&cp, buf, "Gene-commentary", "EC", "EC_number", h) ;
      cp = cp0 ; geneXmlExportXmltag2Acetag (&cp, buf, "Prot-ref_name_E", "Properties", TRUE, h) ;
      cp = cp0 ; geneXmlExportXmltag2Acetag (&cp, buf, "Gene-ref_syn_E", "Locus", TRUE, h) ;
      cp = cp0 ; geneXmlExportXmltag2Acetag (&cp, buf, "Maps_display-str", "Cytogenetic", TRUE, h) ;
      cp = cp0 ; geneXmlExportXmlTagWithValue2Acetag (geneid, &cp, buf, "Gene-commentary", "phenotype", "Locus_Phenotype", h) ;
      cp = cp0 ; generef = xmlGetNextTagContent (&cp, "Gene-ref", h) ;
      if (generef)
	{
	  while ((cp = xmlGetNextTagContent (&generef, "Dbtag_db", h)))
	    {
	      if (!strcmp (cp, "MIM"))
		{
		  cq = xmlGetNextTagContent (&generef, "Object-id_id", h) ;
		  if (cq)
		    {
		      vtxtPrintf (buf, "\nExtern \"OMIM_%s\"\nOMIM\n", cq) ;
		      vtxtPrintf (buf, "\nGeneId %s\nOMIM_mol \"OMIM_%s\" \n", geneid, cq) ;
		    }
		}
	      else if (!strcmp (cp, "MGI"))
		{
		  cq = xmlGetNextTagContent (&generef, "Object-id_id", h) ;
		  if (cq)
		    {
		      vtxtPrintf (buf, "\nExtern \"MGI_%s\"\nMGI\n", cq) ;
		      vtxtPrintf (buf, "GeneId %s\n\n", geneid) ;
		      vtxtPrintf (buf, "GeneId %s\n", geneid) ;
		    }
		}
	    }	  
	}
      cp = cp0 ; geneXmlExportKegg (geneid, &cp, buf, h) ;
      cp = cp0 ; geneXmlExportAccession (&cp, buf, h) ;
    }
  return vtxtPtr (buf) ;
} /* geneXmlGeneContent */

/*************************************************************/

static int geneXmlGene (vTXT vtxt, int line, const char *target)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT buf = vtxtHandleCreate (h) ;
  int ng = 0 ;
  char *cp0, *cp, *geneid ;
  
  cp = cp0 = vtxtPtr (vtxt) ;
  cp = cp0 ; geneid = xmlGetNextTagContent (&cp, "Gene-track_geneid", h) ;

  if (target)
    {
      if (!strcasecmp (target, geneid))
	{
	  printf ("%s\n", cp0) ; 
	  exit (0) ;
	}
    }
  else if (geneid)
    {
      cp = geneXmlGeneContent (geneid, cp0, buf, line, h) ;
      if (cp)
	{ 
	  int n2 = 0 ;
	  messAllocMaxStatus (&n2) ; 
	  ng++; 
	  printf ("geneId %s // %d  maxMem=%d\n%s\n\n", geneid, line,n2, cp) ;
	}
    }
  ac_free (h) ;
  
  return ng ;
} /* geneXmlGene */
     
/*************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  int level = 0, line = 0, oldLine = 0, ng = 0 ;
  const char *target = 0 ;
  vTXT vtxt = vtxtHandleCreate (h) ;
  char *cp ;

  freeinit () ;
  
  getCmdLineOption (&argc, argv, "-t", &target) ;
  if (!target)
    fprintf (stderr, "// Convert stdin to .ace file on stdout, except if  [-t geneid] then export one .xml geneid \n\n") ;
  
  level = freesetfile (stdin, 0) ;
  freespecial ("\n") ;
  while (freecard (level))
    {
      line++ ;
      cp = freepos () ;
      vtxtPrintf (vtxt, "%s\n", cp) ;
      if (strstr (cp, "<Entrezgene-Set>"))
	vtxtClear (vtxt) ;
      if (strstr (cp, "</Entrezgene>"))
	{
	  ng += geneXmlGene (vtxt, oldLine, target) ;
	  oldLine = line ;
	  vtxtClear (vtxt) ;
	}
    }
  fprintf (stderr,"// read %d lines \n", line) ;
      
  ac_free (h) ;
  fprintf (stderr, "// read %d lines, exported %d genes \n", line, ng) ;
  fprintf (stderr, "// done\n") ;
  return 0 ;
}

