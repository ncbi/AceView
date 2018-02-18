#include "misc.h"

#include <sys/types.h>
#include <sys/stat.h>
/* Alliant's Concentrix <sys/stat.h> is hugely deficient */
/* Define things we require in this program              */
/* Methinks S_IFMT and S_IFDIR aren't defined in POSIX   */
#ifndef S_ISDIR
#define S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#endif /*!S_ISDIR*/
#ifndef S_ISREG
#define S_ISREG(m)      (((m)&S_IFMT) == S_IFREG)
#endif /*!S_ISREG*/

int is_directory(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

int is_file(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISREG(buf.st_mode);
}

int file_exists(char * fn)
{
    struct stat buf;
    return ( stat(fn,&buf) == 0);
}

int file_size(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) != 0) return 0;
    return buf.st_size;
}

/*
 * ---------------------------------------------------------------------------
 * File of filename management
 */

FILE *open_fofn(char *files) {
  return fopen(files, "r");
}

char *read_fofn(FILE *fp) {
  char line[256];
  static char name[256];
  
  while (fgets(line, 254, fp)) {
    if (1 == sscanf(line, "%s", name))
      return name;
  }

  return NULL;
}

void close_fofn(FILE *fp) {
  fclose(fp);
}
