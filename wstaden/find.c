#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *myfind(char *file, char* searchpath, int (*found) (char *) )
{
    static char wholePath[1024];
    char *path;
    char *delimiters=":";
    char *f;

    f = NULL;
    if (found(file)) {
	strcpy(wholePath,file);
	f = wholePath;
    } else if (searchpath != NULL) {
	char *paths;

	paths = (char *) malloc(strlen(searchpath)+1);
	strcpy(paths,searchpath);

	path = (char *) strtok(paths,delimiters);
	while (path!= NULL) {

	    (void) strcpy(wholePath,path);
	    (void) strcat(wholePath,"/");
	    (void) strcat(wholePath,file);
	    if (found(wholePath)) {
		f = wholePath;
		break;
	    }
	    path = (char *) strtok((char *)NULL,delimiters);
	}
	free(paths);
    }

    return f;
}
