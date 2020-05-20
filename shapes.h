/*
 * shapes.h
 *      Header file for shapes.c
 *      This file contains methods to create a variety of shapes.
 *          - Circle    (lab03)
 *          - Cone      (lab03)
 *          - Sphere    (Project 3 / Lab07)
 *          - Torus     (Project 1)
 *          - Spring    (Project 1)
 * 
 *  Created on: Sep 17, 2019
 *      Author: Jakob Strobl
 * 
 */
#ifndef _SHAPES_H_
#define _SHAPES_H_

#include <stdio.h>
#include "matrix.h"

typedef struct {
    vec4 A;     /* Point A */
    vec4 B;     /* Point B */
    vec4 C;     /* Point C */
} triangle;
extern vec4 cube_vertices[36];

/**
 * @brief Get a Rectangle based of the parameters
 * 
 * @param pointA - x, y, z position of the bottom-left corner.
 * @param pointB - x, y, z position of the bottom-right corner.
 * @param pointC - x, y, z position of the top-right corner.
 * @param pointD - x, y, z position of the top-left corner.
 * @return vec4* - Vertex array of 2 triangles (6 vertices) to form the rectangle
 */
vec4* getRectangle(vec4 pointA, vec4 pointB, vec4 pointC, vec4 pointD);


/**
 * @brief Create a rectanular prism centered at the origin.
 *  This tries to do it in a weird mathematical way using rotation
 * 
 * @param length - size along the y-axis
 * @param width - size along the x-axis
 * @param depth - size along the z-axis
 * @return vec4* - Vertices to form a prism
 */
vec4* createPrismMath(GLfloat length, GLfloat width, GLfloat depth);

/**
 * @brief Create a rectanular prism centered at the origin.
 *  This tries to do it in a wacky way using pointers
 * 
 * @param length - size along the y-axis
 * @param width - size along the x-axis
 * @param depth - size along the z-axis
 * @return vec4* - Vertices to form a prism
 */
vec4* createPrismWackyPointers(GLfloat length, GLfloat width, GLfloat depth);

/**
 * @brief Create a rectanular prism centered at the origin.
 * 
 * @param length - size along the y-axis
 * @param width - size along the x-axis
 * @param depth - size along the z-axis
 * @return vec4* - Vertices to form a prism
 */
vec4* createPrism(GLfloat length, GLfloat width, GLfloat depth);

vec4* getRubiksCube();

vec4* getSphere(int numRings, int numSides);

// Project 1
vec4* renderTorus(GLfloat R, GLfloat r, int numRings, int numSides);
vec4* renderSpring(GLfloat R, GLfloat r, int numRings, int numSides, int numRotations);

/** @funtion renderCircle
 *      Calculate a circle on 2D plane (x and y axis; z axis = 0) by using the given number of triangles. 
 *      You can also modify the size by changing the radius 
 * 
 *  @param int numTriangles 
 *      The number of triangles we should use to render the circle.
 * 
 *  @return  vec4 pointer (vec4 array) 
 *      Contains the position data of each triangle to form a circle. 
 *      We need 3 vertices per triangle. So the size of vec4 array is:
 *              vec4[numTriangles * 3]
 */ 
vec4* renderCircle(int numTriangles);

/** @funtion renderCone
 *      Calculate a 3D cone, where the base sits on the x-z plane and points up parallel to y axis.
 *      The origin is at the center of the cone.
 *      You can also modify the size of the cone, by changing parameter size.
 * 
 *  @param int numTriangles 
 *      The number of triangles we should use to render the circle.
 * 
 *  @return vec4 * (vec4 array)
 *      Contains the position data of each triangle to form a cone.
 *      We need 3 vertices per triangle, so the size of vec4 array is:
 *          vec4[numTriangles * 3]
 */ 
vec4* renderCone(int numTriangles);

#endif