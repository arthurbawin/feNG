/*
 *  glfem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include "feNG.h"

#if defined(HAVE_GLFW)
#include <../contrib/glfw/include/GLFW/glfw3.h>
#endif

void glfemDrawNodes(double *x, double *y, int n);
void glfemDrawCurve(double *x, double *y, int n);
void glfemDrawCurveDiscrete(double *x, double *y, int n);
void glfemDrawCurveDG(double *x, double *y, int n, int m);
void glfemReshapeWindowsBox(double minX, double maxX, double minY, double maxY, int w, int h);

void glfemMessage(char *message);
void glfemDrawMessage(int h, int v, char *message);
void glfemSetRasterSize(int width, int height);

#if defined(HAVE_GLFW)
GLFWwindow *glfemInit(const char *windowName);
#endif

#endif