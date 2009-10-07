/**
 * Test multi-texturing with GL shading language.
 *
 * Copyright (C) 2008  Brian Paul   All Rights Reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * BRIAN PAUL BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */



#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GL/glut.h"
#include "readtex.h"
#include "extfuncs.h"
#include "shaderutil.h"

static const char *Demo = "multitex";

static const char *VertFile = "multitex.vert";
static const char *FragFile = "multitex.frag";

static const char *TexFiles[2] = 
   {
      "../images/tile.rgb",
      "../images/tree2.rgba"
   };


static GLuint Program;

static GLfloat Xrot = -90.0, Yrot = .0, Zrot = 0.0;
static GLfloat EyeDist = 10;
static GLboolean Anim = GL_TRUE;


/* value[0] = tex unit */
static struct uniform_info Uniforms[] = {
   { "tex1",  1, GL_INT, { 0, 0, 0, 0 }, -1 },
   { "tex2",  1, GL_INT, { 1, 0, 0, 0 }, -1 },
   END_OF_UNIFORMS
};


static void
DrawPolygon(GLfloat size)
{
   glPushMatrix();
   glRotatef(90, 1, 0, 0);
   glNormal3f(0, 0, 1);
   glBegin(GL_POLYGON);

   glMultiTexCoord2f(GL_TEXTURE0, 0, 0);
   glMultiTexCoord2f(GL_TEXTURE1, 0, 0);
   glVertex2f(-size, -size);

   glMultiTexCoord2f(GL_TEXTURE0, 2, 0);
   glMultiTexCoord2f(GL_TEXTURE1, 1, 0);
   glVertex2f( size, -size);

   glMultiTexCoord2f(GL_TEXTURE0, 2, 2);
   glMultiTexCoord2f(GL_TEXTURE1, 1, 1);
   glVertex2f( size,  size);

   glMultiTexCoord2f(GL_TEXTURE0, 0, 2);
   glMultiTexCoord2f(GL_TEXTURE1, 0, 1);
   glVertex2f(-size,  size);

   glEnd();
   glPopMatrix();
}


static void
draw(void)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glPushMatrix(); /* modelview matrix */
      glTranslatef(0.0, 0.0, -EyeDist);
      glRotatef(Zrot, 0, 0, 1);
      glRotatef(Yrot, 0, 1, 0);
      glRotatef(Xrot, 1, 0, 0);

      DrawPolygon(3.0);

   glPopMatrix();

   glutSwapBuffers();
}


static void
idle(void)
{
   GLfloat t = 0.05 * glutGet(GLUT_ELAPSED_TIME);
   Yrot = t;
   glutPostRedisplay();
}


static void
key(unsigned char k, int x, int y)
{
   (void) x;
   (void) y;
   switch (k) {
   case ' ':
   case 'a':
      Anim = !Anim;
      if (Anim)
         glutIdleFunc(idle);
      else
         glutIdleFunc(NULL);
      break;
   case 'z':
      EyeDist -= 0.5;
      if (EyeDist < 3.0)
         EyeDist = 3.0;
      break;
   case 'Z':
      EyeDist += 0.5;
      if (EyeDist > 90.0)
         EyeDist = 90;
      break;
   case 27:
      exit(0);
   }
   glutPostRedisplay();
}


static void
specialkey(int key, int x, int y)
{
   GLfloat step = 2.0;
   (void) x;
   (void) y;
   switch (key) {
   case GLUT_KEY_UP:
      Xrot += step;
      break;
   case GLUT_KEY_DOWN:
      Xrot -= step;
      break;
   case GLUT_KEY_LEFT:
      Yrot -= step;
      break;
   case GLUT_KEY_RIGHT:
      Yrot += step;
      break;
   }
   glutPostRedisplay();
}


/* new window size or exposure */
static void
Reshape(int width, int height)
{
   GLfloat ar = (float) width / (float) height;
   glViewport(0, 0, (GLint)width, (GLint)height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-2.0*ar, 2.0*ar, -2.0, 2.0, 4.0, 100.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}


static void
InitTextures(void)
{
   GLenum filter = GL_LINEAR;
   int i;

   for (i = 0; i < 2; i++) {
      GLint imgWidth, imgHeight;
      GLenum imgFormat;
      GLubyte *image = NULL;

      image = LoadRGBImage(TexFiles[i], &imgWidth, &imgHeight, &imgFormat);
      if (!image) {
         printf("Couldn't read %s\n", TexFiles[i]);
         exit(0);
      }

      glActiveTexture(GL_TEXTURE0 + i);
      glBindTexture(GL_TEXTURE_2D, 42 + i);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 4, imgWidth, imgHeight,
                        imgFormat, GL_UNSIGNED_BYTE, image);
      free(image);
      
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
   }
}


static GLuint
CreateProgram(const char *vertProgFile, const char *fragProgFile,
              struct uniform_info *uniforms)
{
   GLuint fragShader, vertShader, program;

   vertShader = CompileShaderFile(GL_VERTEX_SHADER, vertProgFile);
   fragShader = CompileShaderFile(GL_FRAGMENT_SHADER, fragProgFile);
   assert(vertShader);
   program = LinkShaders(vertShader, fragShader);

   glUseProgram_func(program);

   InitUniforms(program, uniforms);

   return program;
}


static void
InitPrograms(void)
{
   Program = CreateProgram(VertFile, FragFile, Uniforms);
}


static void
InitGL(void)
{
   const char *version = (const char *) glGetString(GL_VERSION);

   if (version[0] != '2' || version[1] != '.') {
      printf("Warning: this program expects OpenGL 2.0\n");
      /*exit(1);*/
   }
   printf("GL_RENDERER = %s\n",(const char *) glGetString(GL_RENDERER));

   GetExtensionFuncs();

   InitTextures();
   InitPrograms();

   glEnable(GL_DEPTH_TEST);

   glClearColor(.6, .6, .9, 0);
   glColor3f(1.0, 1.0, 1.0);
}


int
main(int argc, char *argv[])
{
   glutInit(&argc, argv);
   glutInitWindowSize(500, 400);
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutCreateWindow(Demo);
   glutReshapeFunc(Reshape);
   glutKeyboardFunc(key);
   glutSpecialFunc(specialkey);
   glutDisplayFunc(draw);
   if (Anim)
      glutIdleFunc(idle);
   InitGL();
   glutMainLoop();
   return 0;
}
