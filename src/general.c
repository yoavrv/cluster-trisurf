/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdio.h>
#include<stdlib.h>
#include "general.h"
#include<stdarg.h>

#include <sys/time.h>
#include <unistd.h>
#include <time.h>

#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>

ts_uint ts_fprintf(FILE *fd, char *fmt, ...){
if(quiet) return TS_SUCCESS;
	va_list ap;
	va_start(ap,fmt);
	char tmbuf[255];
	struct timeval now;
  	gettimeofday(&now, 0);
	strftime(tmbuf, sizeof tmbuf, "%Y-%m-%d %H:%M:%S", localtime(&now.tv_sec));
fprintf(fd, "[%s] ",tmbuf); 
vfprintf(fd, fmt, ap); /* Call vfprintf */
va_end(ap); /* Cleanup the va_list */
return TS_SUCCESS;
}

void err(char *text){
	ts_fprintf(stderr,"Err: %s\n", text);
}

void fatal(char *text, ts_int errcode){
	ts_fprintf(stderr,"Fatal: %s. TERMINATED!\n", text);
	exit(errcode);
}


/* Open/create the file named in 'pidFile', lock it, optionally set the
   close-on-exec flag for the file descriptor, write our PID into the file,
   and (in case the caller is interested) return the file descriptor
   referring to the locked file. The caller is responsible for deleting
   'pidFile' file (just) before process termination. 'progName' should be the
   name of the calling program (i.e., argv[0] or similar), and is used only for
   diagnostic messages. If we can't open 'pidFile', or we encounter some other
   error, then we print an appropriate diagnostic and terminate. */

/* 
 This is filelock/create_pid_file.c (Listing 55-4, page 1143), an example program file from the book, The Linux Programming Interface.

The source code file is copyright 2010, Michael Kerrisk, and is licensed under the GNU Lesser General Public License, version 3. 
*/

#define BUF_SIZE 100
int createPidFile(const char *progName, const char *pidFile, int flags)
{
    int fd;
    char buf[BUF_SIZE];

    fd = open(pidFile, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    if (fd == -1){
        ts_fprintf(stderr,"Could not open PID file %s", pidFile);
	fatal("Cannot continue (1)",1);
}
    if (flags & CPF_CLOEXEC) {

        /* Set the close-on-exec file descriptor flag */

        /* Instead of the following steps, we could (on Linux) have opened the
           file with O_CLOEXEC flag. However, not all systems support open()
           O_CLOEXEC (which was standardized only in SUSv4), so instead we use
           fcntl() to set the close-on-exec flag after opening the file */

        flags = fcntl(fd, F_GETFD);                     /* Fetch flags */
        if (flags == -1){
            ts_fprintf(stderr,"Could not get flags for PID file %s", pidFile);
	fatal("Cannot continue (2)",1);
}
        flags |= FD_CLOEXEC;                            /* Turn on FD_CLOEXEC */

        if (fcntl(fd, F_SETFD, flags) == -1)            /* Update flags */
            ts_fprintf(stderr,"Could not set flags for PID file %s", pidFile);
		fatal("Cannot continue (3)",1);
	    
    }

    if (lockRegion(fd, F_WRLCK, SEEK_SET, 0, 0) == -1) {
        if (errno  == EAGAIN || errno == EACCES){
            ts_fprintf(stderr,"PID file '%s' is locked; probably "
                     "'%s' is already running", pidFile, progName);
		fatal("Cannot continue (4)",1);
}
        else{
            ts_fprintf(stderr,"Unable to lock PID file '%s'", pidFile);
	fatal("Cannot continue (5)",1);
}
    }

    if (ftruncate(fd, 0) == -1){
        ts_fprintf(stderr,"Could not truncate PID file '%s'", pidFile);
	fatal("Cannot continue (6)",1);
}

    snprintf(buf, BUF_SIZE, "%ld\n", (long) getpid());
    if (write(fd, buf, strlen(buf)) != strlen(buf)){

        ts_fprintf(stderr,"Writing to PID file '%s'", pidFile);
	fatal("Cannot continue (7)",1);
}
    return fd;
}



/* Lock a file region (private; public interfaces below) */

static int
lockReg(int fd, int cmd, int type, int whence, int start, off_t len)
{
    struct flock fl;

    fl.l_type = type;
    fl.l_whence = whence;
    fl.l_start = start;
    fl.l_len = len;

    return fcntl(fd, cmd, &fl);
}


int                     /* Lock a file region using nonblocking F_SETLK */
lockRegion(int fd, int type, int whence, int start, int len)
{
    return lockReg(fd, F_SETLK, type, whence, start, len);
}

