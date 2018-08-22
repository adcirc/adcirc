#include "cfi.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

void f_makedir(char* path, int len)
{
  char *dirnm;

  dirnm = (char *) malloc(len+1);
  memcpy(dirnm, path, len);

  dirnm[len] = '\0';

#ifdef _WIN32
  CreateDirectory(dirnm,NULL);
#else
  mkdir(dirnm, 0755);
#endif


  free(dirnm);
}

void MAKEDIR(protoFSTRING(fpath) protoLenFSTRING(fpath))
{
  char *path  = FCD2CP(fpath);
  int   len   = FCDLEN(fpath);

  f_makedir(path, len);
}

void makedir(protoFSTRING(fpath) protoLenFSTRING(fpath))
{
  char *path  = FCD2CP(fpath);
  int   len   = FCDLEN(fpath);

  f_makedir(path, len);
}

void makedir_(protoFSTRING(fpath) protoLenFSTRING(fpath))
{
  char *path  = FCD2CP(fpath);
  int   len   = FCDLEN(fpath);

  f_makedir(path, len);
}
