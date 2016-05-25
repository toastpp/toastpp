#ifndef EGLCONFIG_INCLUDED
#define EGLCONFIG_INCLUDED


#include "egltypedefs.h"
#include <GLES/gl.h>


#define MAX_ATTRIBS 128
#define FIRST_ATTRIB EGL_BUFFER_SIZE


struct _egl_config
{
   EGLConfig Handle;   /* the public/opaque handle which names this config */
   EGLint Attrib[MAX_ATTRIBS];
};


#define SET_CONFIG_ATTRIB(CONF, ATTR, VAL) \
   assert((ATTR) - FIRST_ATTRIB < MAX_ATTRIBS); \
   ((CONF)->Attrib[(ATTR) - FIRST_ATTRIB] = VAL)


#define GET_CONFIG_ATTRIB(CONF, ATTR) ((CONF)->Attrib[(ATTR) - FIRST_ATTRIB])


extern void
_eglInitConfig(_EGLConfig *config, EGLint id);


extern EGLConfig
_eglGetConfigHandle(_EGLConfig *config);


extern _EGLConfig *
_eglLookupConfig(_EGLDriver *drv, EGLDisplay dpy, EGLConfig config);


extern _EGLConfig *
_eglAddConfig(_EGLDisplay *display, _EGLConfig *config);


extern EGLBoolean
_eglParseConfigAttribs(_EGLConfig *config, const EGLint *attrib_list);


extern EGLBoolean
_eglChooseConfig(_EGLDriver *drv, EGLDisplay dpy, const EGLint *attrib_list, EGLConfig *configs, EGLint config_size, EGLint *num_config);


extern EGLBoolean
_eglGetConfigAttrib(_EGLDriver *drv, EGLDisplay dpy, EGLConfig config, EGLint attribute, EGLint *value);


extern EGLBoolean
_eglGetConfigs(_EGLDriver *drv, EGLDisplay dpy, EGLConfig *configs, EGLint config_size, EGLint *num_config);


extern void
_eglSetConfigAttrib(_EGLConfig *config, EGLint attr, EGLint val);


#endif /* EGLCONFIG_INCLUDED */
