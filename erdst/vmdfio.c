/* ERmod - Energy Representation Module
   Copyright (C) 2000- The ERmod authors

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
   As a special exception, you may use this file as part of a free software
   without restriction.  Specifically, if other files instantiate
   templates or use macros or inline functions from this file, or you compile
   this file and link it with other files to produce an executable, this
   file does not by itself cause the resulting executable to be covered by
   the GNU General Public License.  

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. */

/*
  Generic trajectory access via VMD plugin.
 */

/* FIXME: add F77_FUNC(small,CAPITAL) issues */
#include "config.h"

#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <math.h>

#include <dlfcn.h>

#include "../vmdplugins/molfile_plugin.h"

/*
  Declaration of plugin access points.
*/

#ifndef INSTALL_PLUGIN_PATH
#warning INSTALL_PLUGIN_PATH is not defined
#define INSTALL_PLUGIN_PATH ""
#endif

#define ADHOC_RESOLVER "adhocresolve.list"

#define ERMOD_FORCE_PLUGIN "ERMOD_FORCE_PLUGIN_TYPE"

/* lazy */
#define MAXTYPES 1000
static int typecounts = 0;
static molfile_plugin_t* vmdpluginentries[MAXTYPES];

/* lazy too */
#define MAXPLUGINS 100
static int plugincounts = 0;
static void* ldhandles[MAXPLUGINS];

/*
  Callback function. Must be compatible with vmdplugin_register_cb. 
 */
static int register_callback(void* __vp, vmdplugin_t * entity)
{
  if(entity == NULL || entity -> type == NULL) return 0;
  vmdpluginentries[typecounts++] = (molfile_plugin_t*)entity;
  return 0;
}

static int check_directory_available(char *path)
{
  struct stat st;
  int r = stat(path, &st);

  if(r == 0 && S_ISDIR(st.st_mode))
    return 1;
  return 0;
}

static int is_executable(char* path)
{
  struct stat st;
  int r;
  r = stat(path, &st);
  if(r != 0) return 0;
  if(!(S_ISREG(st.st_mode) || S_ISLNK(st.st_mode))) {
    return 0;
  }
  
  if(access(path, X_OK) == 0) return 1;
  return 0;
}

static char* strdupcat(char* str1, char* str2)
{
  char* buf = malloc(strlen(str1) + strlen(str2) + 1);
  strcpy(buf, str1);
  strcat(buf, str2);
  return buf;
}

static char* find_plugin_paths()
{
  char *pluginpath;
  /* contract: return malloc'ed string, or return NULL if error */

  /* is ERMOD_PLUGINS present and directory? */
  if((pluginpath = getenv("ERMOD_PLUGINS")) && 
     check_directory_available(pluginpath)) {
    return strdup(pluginpath);
  }

  /* is INSTALL_PLUGIN_PATH available? */
  if(strlen(INSTALL_PLUGIN_PATH) > 0 &&
     check_directory_available(INSTALL_PLUGIN_PATH)) {
    return strdup(INSTALL_PLUGIN_PATH);
  }

  return NULL;
}


void vmdfio_init_traj_(void)
{
  char* plugindir;
  
  if((plugindir = find_plugin_paths()) == NULL){
    fprintf(stderr, "Error: cannot find plugin paths, please set ERMOD_PLUGINS environment variable\n");
    exit(EXIT_FAILURE);
  }

  /* traverse directory */
  do{
    DIR* dirp;
    struct dirent *ent;
    dirp = opendir(plugindir);
    
    while((ent = readdir(dirp)) != NULL){
      /* concat paths */
      char* d1 = strdupcat(plugindir, "/");
      char* pluginpath = strdupcat(d1, ent -> d_name);
      free(d1);
      
      if(is_executable(pluginpath) && pluginpath[strlen(pluginpath)-1] != '.'){
	void* handle;
	handle = dlopen(pluginpath, RTLD_NOW | RTLD_GLOBAL);
	if(!handle){
	  fprintf(stderr, "Warning: failed to load plugin \"%s\", reason: \"%s\"\n", pluginpath, dlerror());
	  continue;
	}
	void* initptr;
	void* entptr;
	initptr = dlsym(handle, "vmdplugin_init");
	entptr = dlsym(handle, "vmdplugin_register");
	if(initptr != NULL && entptr != NULL){
	  int r;
	  r = ((int (*)(void))initptr)();
	  if(r != 0)
	    fprintf(stderr, "Warning: error while initializing %s\n", pluginpath);
	  ldhandles[plugincounts++] = handle;
	  r = ((int (*)(void*, vmdplugin_register_cb))entptr)(NULL, register_callback);
	  if(r != 0)
	    fprintf(stderr, "Warning: error while registering %s\n", pluginpath);
	}else{
	  fprintf(stderr, "Warning: failed to load entry point (plugin: %s)\n", pluginpath);
	}
      }
    }
  }while(0);
}

void vmdfio_fini_traj_(void)
{
  int i;
  for(i = 0; i < plugincounts; ++i){
    void* entptr = dlsym(ldhandles[i], "vmdplugin_fini");
    if(entptr){
      ((int (*)(void))entptr)();
    }
    dlclose(ldhandles[i]);
  }
}

typedef struct vmdpluginio_t {
  molfile_plugin_t * plugin;
  void* filehandle;
  int natoms;
} vmdpluginio;

void vmdfio_open_traj_(void **handle, char *fname, int *fnamelen, int *status)
{
  char* buf;
  size_t buflen = 8192;
  ssize_t r;
  char* ext;
  int i;

  buf = malloc(sizeof(char) * buflen + 1);
  strncpy(buf, fname, *fnamelen);
  buf[*fnamelen] = '\0';

  while(1){
    r = readlink(buf, buf, buflen);
    if(r == -1) break;
    buf[r] = '\0';
  }

  if(*fnamelen == strlen(buf) && strncmp(fname, buf, *fnamelen) == 0){
    /* not a symbolic link? */
    fprintf(stderr, "Error: vmdfio.c: failed to open with vmdfio_open_traj_. (filename = \"%s\".) Perhaps it's not a symbolic link?\n", buf);
    *status = -1;
    goto cleanup;
  }
  
  /* 
     select plugin to use.
     use filename_extension to select
  */
  *status = -1;
  *handle = NULL;
  {
    /* get filename extension */
    char* lastdot = strrchr(buf, '.');

    if(lastdot == NULL) goto cleanup;
    ext = lastdot + 1;
  }
  fprintf(stderr, "Opening: \"%s\"...\n", buf);

  if(getenv(ERMOD_FORCE_PLUGIN) != NULL && 
     strlen(getenv(ERMOD_FORCE_PLUGIN)) > 0) {
    ext = getenv(ERMOD_FORCE_PLUGIN);
  }

  for(i = 0; i < typecounts; ++i){
    int extlen = strlen(ext);
    molfile_plugin_t *p = vmdpluginentries[i];
    char* plugin_supportext;
    char* tokptr;
    char* ptr;

    if(strcmp(p -> type, MOLFILE_PLUGIN_TYPE) != 0) continue; 
    plugin_supportext = strdup(p -> filename_extension);

    while(1){
      ptr = strtok_r(plugin_supportext, ",", &tokptr);
      if(ptr == NULL) break;
      plugin_supportext = NULL; /* in the next subsequent call NULL is passed */

      if(strncmp(ptr, ext, extlen) == 0 && 
	 (ptr[extlen] == '\0' || ptr[extlen] == ',')){
	void* fh;
	vmdpluginio* pp;

	fprintf(stderr, "  Trying plugin \"%s\"...", p -> prettyname);

	if(p -> abiversion <= 10) {
	  fprintf(stderr, "Error: \"%s\" supports trajectory format, but it is too old (requires ABI version > 10)\n", p -> prettyname);
	  fprintf(stderr, "Please update plugin files to those of the latest VMD's\n");
	  continue;
	}

	pp = malloc(sizeof(vmdpluginio));
	pp -> plugin = p;

	/* Found, open with this plugin */
	pp -> natoms = MOLFILE_NUMATOMS_UNKNOWN;
	fh = (p -> open_file_read)(buf, ext, &(pp->natoms));
	if(fh == NULL){
	  fprintf(stderr,
		  "Error!\n"
		  "Plugin supports file format \"%s\", estimated from the extension \"%s\".\n"
		  "But the file format of \"%s\" does not match to the estimated format.\n"
		  "This is possibly that the file has an incorrect extension, the file is a nested symbolic link (which is unsupported), or the trajectory file header is corrupt.\n",
		  p -> prettyname, ext, buf);
	  /* Special error message for AMBER */
	  if(strcmp(ext, "nc") == 0) {
	    fprintf(stderr,
		    "Perhaps you forgot to specify \"ioutfm = 1\" in AMBER?\n");
	  }
	  exit(1);
	}
	pp -> filehandle = fh;
	fprintf(stderr, "OK\n");

	*handle = pp;
	break;
      }
    }
    free(plugin_supportext);
    if(*handle){
      *status = 0;
      break;
    }
  }

  if(*status == -1){
    fprintf(stderr, "Error: could not find appropriate plugins\n");
  }


 cleanup:
  free(buf);
  return;
}

void vmdfio_read_traj_step_(void **handle, float* xout, float* box, int *natoms_aux, int *status)
{
  vmdpluginio *p = *handle;
  molfile_plugin_t *plugin = p -> plugin;
  molfile_timestep_t snapshot;
  int natoms = p -> natoms;
  float* buf;
  int r;

  if(natoms == MOLFILE_NUMATOMS_UNKNOWN){
    /* not determined from the trajectory */
    natoms = *natoms_aux;
  }else{
    /* check integrity */
    if(natoms != *natoms_aux){
      fprintf(stderr, "Error: # of atoms in trajectory does not match with # of atoms in configurations. Perhaps you mistook the trajectory?\n");
      *status = -1;
      return;
    }
  }

  buf = malloc(sizeof(float) * 3 * natoms);
  snapshot.coords = buf;

  r = (plugin -> read_next_timestep)(p -> filehandle, natoms, &snapshot);

  do{
    int i;
    double x, y, u, v, w;

    if(r == MOLFILE_EOF){
      *status = -1;
      break;
    }

    /* both are the same format */
    for(i = 0; i < natoms * 3; ++i)
      xout[i] = buf[i];

    /*
      ~a = (1, 0, 0)
      ~b = (x, y, 0)
      ~c = (u, v, w)
      ~a.~b = x = cos gamma
      |~a*~b| = y = sin gamma
      ~a.~c = u = cos beta
      ~b.~c = xu + yv = cos alpha
    */

    x = cos(snapshot.gamma * M_PI / 180.0);
    y = sin(snapshot.gamma * M_PI / 180.0);
    u = cos(snapshot.beta * M_PI / 180.0);
    v = (cos(snapshot.alpha * M_PI / 180.0) - x * u) / y; /* FIXME: potential underflow risk */
    w = sqrt(1 - u * u - v * v);  /* FIXME: same above */

    box[0] = snapshot.A;
    box[3] = snapshot.B * x; 
    box[4] = snapshot.B * y;
    box[6] = snapshot.C * u; 
    box[7] = snapshot.C * v;
    box[8] = snapshot.C * w;
    
    box[1] = 0.; box[2] = 0.; box[5] = 0;

    *status = 0;
  }while(0);

  free(buf);
  return;
}

void vmdfio_close_traj_(void **handle)
{
  vmdpluginio *p = *handle;
  molfile_plugin_t *plugin = p -> plugin;
  (plugin -> close_file_read)(p -> filehandle);
  free(*handle);
}
