#include <wac/ac.h>

int main(int argc, char * argv[],char * envp[])
{   
  if (! argv[1])
    {
      fprintf(stderr,"must give -db or -treat\n");
      exit(1);
    }
  
  if (strcmp(argv[1],"-db") == 0)
    {
      AC_DB db;
      int num_authors=1;
      char *query ; const char *s;
      AC_ITER it;
      AC_OBJ ob;
      
      query = " find paper ?cgc* ";
      if ( argv[3] )
	query = argv[3];
      
      db = ac_open_db (argv[2],&s);
      if (!db)
	{
	  fprintf(stderr,"can't open database: %s\n",s);
	  exit(1);
	}
      
      it = ac_query_iter (db, 0, query, NULL, NULL);
      
      while ((ob = ac_next_obj (it)))
	{
	  const char *j, *p, *v;
	  int y;
	  AC_TABLE authors;
	  AC_OBJ tobj;
	  
	  tobj = ac_tag_obj(ob,"Journal",NULL);
	  if (tobj)
	    {
	      j = ac_name(tobj);
	      ac_free(tobj);
	    }
	  else
	    j = 0;
	  if (!j || (j && (strstr(j,"Worm") || strstr(j,"Meeting"))))	
	    continue;
	  
	  y=ac_tag_int(ob,"Year",0);
	  
	  v=ac_tag_text(ob,"Volume",0);
	  
	  p=ac_tag_text(ob,"Page",0);
	  
	  authors=ac_tag_table(ob,"Author",NULL);
	  
	  if (! authors)
	    {
	      fprintf(stderr,"%s - no authors!\n",ac_name(ob));
	    }
	  
	  if(!j && !y && !v && !p && !authors)continue;
	  if(j)printf("%s|" ,j);else printf("|");
	  if(y)printf("%d|" ,y);else printf("|");
	  if(v)printf("%s|" ,v);else printf("|");
	  if(p)printf("%s|" ,p);else printf("|");
	  if(authors)
	    {
	      int ir;
	      for (ir=0; ir < num_authors && ir < authors->rows; ir++)
		{
		  tobj = ac_table_obj(authors, ir, 0, NULL);
		  if (!tobj)
		    break;
		  s = ac_name(tobj);
		  printf("%s%s",ir? "," : "", s);
		  ac_free(tobj);
		}
	    } 
	  printf("|%s|\n",ac_name(ob));
	  
	  ac_free(authors);
	  ac_free(ob);
	  ac_db_close (db) ;
	}
    }
  else if (strcmp(argv[1],"-treat") == 0)
    {
      char b[1000];
      
      while (fgets(b,sizeof(b),stdin))
	{
	  char *pmid, *id, *s;
	  pmid = strrchr(b,'|');
	  if (!pmid)
	    continue;
	  *pmid++ = 0;
	  s = strchr(pmid,'\n');
	  if (s) 
	    *s = 0;
	  id = strrchr(b,'|');
	  if (!id)
	    continue;
	  *id++ = 0;
	  if (atoi(pmid))
	    {
	      printf("Paper : pm%s\n",pmid);
	      printf("CGCID \"%s\"\n",id);
	    }
	  else
	    {
	      printf("// Paper : \"???\"  // %s\n",pmid);
	      printf("// CGCID \"%s\"\n",id);
	    }
	  printf("\n");
	}
    }
  
  return 0;
}

