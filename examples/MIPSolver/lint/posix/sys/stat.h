
/*  sys/stat.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__STAT_H)
#define __STAT_H

struct  stat {
	dev_t   st_dev;
	ino_t   st_ino;
	mode_t  st_mode;
	nlink_t st_nlink;
	uid_t   st_uid;
	gid_t   st_gid;
	off_t   st_size;
	time_t  st_atime;
	time_t  st_mtime;
	time_t  st_ctime;
};

#define  S_ISDIR  0040000
#define  S_ISCHR  0020000
#define  S_ISBLK  0060000
#define  S_ISREG  0100000
#define  S_ISFIFO 0010000

#define S_ISUID 04000
#define S_ISGID 02000
#define S_IRWXU 00700
#define S_IRUSR 00400
#define S_IWUSR 00200
#define S_IXUSR 00100
#define S_IRWXG 00070
#define S_IRGRP 00040
#define S_IWGRP 00020
#define S_IXGRP 00010
#define S_IRWXO 00007
#define S_IROTH 00004
#define S_IWOTH 00002
#define S_IXOTH 00001

int chmod(const char *, mode_t );
int fstat(int, struct stat *);
int mkdir(const char *, mode_t );
int mkfifo(const char *, mode_t );
int stat(const char *, struct stat * );
mode_t umask(mode_t );

#endif
