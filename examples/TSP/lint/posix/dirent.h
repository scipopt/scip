
/*  dirent.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__DIRENT_H)
#define __DIRENT_H

struct dirent
	{
	char        d_name[1];
};



typedef struct
	{
	int _DIR1;
}       DIR;


DIR           *opendir(const char * );
struct dirent *readdir(DIR * );
void           rewinddir(DIR * );
int            closedir(DIR * );


#endif
