#include    <wac/ac.h>

/*****************************************************************************************/
/* This code is part of the AceView Genome annotation Pipeline 
 * Author mieg@ncbi, may 2003
 * 
 * The problem is that in locuslink, the GO annotation is given by name
 * and number, but not separated in biol, molec and cellular
 * here we read the data in LocusLinkInfo
 * in the tag Gene_ontology Text Tex
 * expecting by example :  Gene_ontology "GO:0030282"      "bone mineralization"
 * and tun it into Go_c_ace since this one is cellular
 */

int main (int argc, char *argv[])
{   
  AC_DB    db ;
  AC_ITER  iter ;
  AC_OBJ   obj, oNm ;
  AC_TABLE tt, tt1 ;
  char *which, *goname ;
  const char *nm, *type, *error ;
  int ir ;
  
  if(argc<2) messcrash ("Usage: \n programName AceDB_path [thename]\n");
  which = (argc >= 3) ? argv[2] : "*" ;
    
  aceInit (argv[1]) ;  /* initialise tace, since we run as a standalone code */
  db=ac_open_db (argv[1], &error) ;
  if(db &&
     (iter=ac_query_iter (db, TRUE, messprintf("find locuslink Gene_ontology AND IS \"%s\"", which), 0, 0)))
    { 
      printf("\n// locuslinkgo extracted\n\n ") ;
      
      while((obj=ac_next_obj (iter)))
	{
	  printf ("LocusLink \"%s\"\n-D Go_b_ace\n-D Go_m_ace\n-D Go_c_ace\n",ac_name(obj) ) ;
	  tt = ac_tag_table (obj, "Gene_ontology", 0) ;
	  
	  for (ir = 0 ; ir < tt->rows ; ir++)
	    {
	      nm = ac_table_printable (tt, ir, 0, 0) ;
	      if (!nm || strncasecmp (nm, "GO:", 3))
		continue ;
	      if ((oNm = ac_get_obj (db, "GO", nm, 0)))
		{
		  if ((tt1 = ac_tag_table (oNm, "Type", 0)))
		    {
		      type = ac_table_printable (tt1, 0, 0, 0) ;
		      goname = strnew (ac_table_printable (tt1, 0, 1, 0), 0) ;
		    
		      if (type && goname)
			{
			  if (!strcasecmp (type, "Molecular_function"))
			    printf ("Go_m_ace \"%s\"\n", goname) ;
			  else if (!strcasecmp (type, "Biological_process"))
			    printf ("Go_b_ace \"%s\"\n", goname) ;
			  else if (!strcasecmp (type, "Cellular_component"))
			    printf ("Go_c_ace \"%s\"\n", goname) ;
			}
		      messfree (goname) ;
		    }
		  ac_free (tt1) ;
		}
	      ac_free (oNm) ;
	    }
	  ac_free (tt) ;
	  printf ("\n\n") ;
	}
    }  
  ac_db_close (db);
  aceQuit (FALSE) ;
  return 0;
}

/*****************************************************************************************/
/*****************************************************************************************/
