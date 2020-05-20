/**
 * @file camera.h
 * @author Jakob Strobl
 * @brief Functions for handling and setting cameras for use in OpenGL
 * @version 0.1
 * @date 2019-10-27
 * 
 * @copyright Copyright (c) 2019
 */
#ifndef CAMERA_H_
#define CAMERA_H_

#include "matrix.h"

typedef struct {
    GLfloat left;
    GLfloat right;
    GLfloat bottom;
    GLfloat top;
    GLfloat near;
    GLfloat far;
} viewVolume; 

/**
 * @brief 
 * 
 * @param eye 
 * @param at 
 * @param up 
 * @return mat4* the model view matrix
 */
mat4* look_at(vec3 eye, vec3 at, vec3 up);

/**
 * @brief 
 * 
 * @param view 
 * @return mat4* - the projection matrix
 */
mat4* ortho(viewVolume* view);

mat4* frustrum(viewVolume* view);

/**
 * @brief 
 * 
 * @param fovy 
 * @param aspect 
 * @param near 
 * @param far 
 * @return mat4* - the perspective matrix
 */
mat4* perspective(GLfloat fovy, GLfloat aspect, GLfloat near, GLfloat far);






#endif