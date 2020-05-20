/**
 * @file camera.c
 * @author Jakob Strobl
 * @brief Functions for handling and setting cameras for use in OpenGL
 * @version 0.1
 * @date 2019-10-27
 * 
 * @copyright Copyright (c) 2019
 */
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "camera.h"

mat4* look_at(vec4 eye, vec4 at, vec4 up) {
    mat4 *model_view = calloc(1, sizeof(mat4));
    vec4 u, v, n; // eye = p
    
    // Calculate n
    //   N is the vector of the eye to what we are looking at. (eye -> at)
    vecSub_inplace(eye, at, n);
    vecNormalize_inplace(n, n); 
    //vecPrint(n);

    // Calculate u
    //   u = norm(v x n);
    vecCrossProd_inplace(up, n, u);
    vecNormalize_inplace(u, u);
    //vecPrint(u);


    // Calculate v
    //   v = norm(n x u);
    vecCrossProd_inplace(n, u, v);
    vecNormalize_inplace(v, v);
    //vecPrint(v);


    model_view->x[0] = u[0];
    model_view->y[0] = u[1];
    model_view->z[0] = u[2];
    
    model_view->x[1] = v[0];
    model_view->y[1] = v[1];
    model_view->z[1] = v[2];

    model_view->x[2] = n[0];
    model_view->y[2] = n[1];
    model_view->z[2] = n[2];

    model_view->w[3] = 1;

    mat4 *rt = matMult(model_view, translate(-eye[0], -eye[1], -eye[2]));
    free(model_view);
    
    return rt;
}

mat4* ortho(viewVolume *view) {
    mat4 *orthographic;

    // Calculate center of viewVolume
    GLfloat x = (view->right + view->left) / 2;
    GLfloat y = (view->top + view->bottom) / 2;
    GLfloat z = (view->near + view->far) / 2;
    
    // Calculate the scaling factor for each axis
    GLfloat x_s = 2 / (view->right - view->left);
    GLfloat y_s = 2 / (view->top - view->bottom);
    GLfloat z_s = 2 / (view->near - view->far);

    // Translate center of mass of the view volume to OpenGL canonical view AND
    //   scale the view volume to OpenGL's canonical view volume.
    return chainCTM(2, translate(-x, -y, -z), scale(x_s, y_s, z_s));
}

mat4* frustrum(viewVolume *view) {
    mat4 *frustrum = calloc(1, sizeof(mat4));

    // Perspective transformation 
    GLfloat x_scale = (-2 * view->near)/(view->right - view->left);
    GLfloat y_scale = (-2 * view->near)/(view->top - view->bottom);

    // Perspective shearing
    GLfloat x_shear = (view->left + view->right)/(view->right - view->left);
    GLfloat y_shear = (view->bottom + view->top)/(view->top - view->bottom);

    // Perspective depth
    GLfloat depth = (view->near + view->far)/(view->far - view->near);
    GLfloat z_offset = ((-2 * view->near * view->far) / (view->far - view->near));

    frustrum->x[0] = x_scale;
    frustrum->y[1] = y_scale;
    frustrum->z[0] = x_shear;
    frustrum->z[1] = y_shear;
    frustrum->z[2] = depth;
    frustrum->z[3] = -1; // Depth?
    frustrum->w[2] = z_offset;

    return frustrum;
}

mat4* perspective(GLfloat fovy, GLfloat aspect, GLfloat near, GLfloat far) {
    GLfloat top = -1 * near * tanf(((fovy * M_PI)/ 180)/2);
    GLfloat right = top/aspect;
    viewVolume view = {-right, right, -top, top, near, far};

    // printf("Top: %f", top);
    // printf("Right: %f", right);

    return frustrum(&view);
}