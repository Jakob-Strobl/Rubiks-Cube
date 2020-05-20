/*
 *  matrix.c
 *      Library for 4x1 Vectors and 4x4 Matrices. 
 *  
 *  @author: Jakob Strobl 
 *  Created on: August 31, 2019 
 */

// Custom definitions based on OS 
#ifdef __APPLE__  // include Mac OS X verions of headers
#define GL_SILENCE_DEPRECATION
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else // non-Mac OS X operating systems
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#endif  // __APPLE__

// Non-OS dependent dependencies 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "matrix.h"

/***      Vector Operations      ***/

/*  @function vecPrint:
 *      Print each element of the vector with 3 points of floating
 *      point granularity
 */
void vecPrint(vec4 vec) {
    printf("[%8.3f, %8.3f, %8.3f, %8.3f]\n", vec[0], vec[1], vec[2], vec[3]);
}

/*  @function vecPrintRaw:
 *      Print each element of the vector with full floating
 *      point granularity (prints complete floating point value).
 */
void vecPrintAll(vec4 vec) {
    printf("[%10f, %10f, %10f, %10f]\n", vec[0], vec[1], vec[2], vec[3]);
}

/*  @function vecScalarMult:
 *      Calculate the Scalar multiplacation of the vector given, and return malloced vec4
 *      Remember to free() the vector when finished. 
 */
vec4* vecScalarMult(GLfloat s, vec4 v) {
    // Allocate space for vec4 type
    vec4* result = malloc(sizeof(vec4));
    GLfloat* ptr = result[0];   // pointer to jump through array 

    *ptr = v[0] * s;
    ptr += 1;                   // Go to index result[1]
    *ptr = v[1] * s;
    ptr += 1;                   // Go to index result[2]
    *ptr = v[2] * s;
    ptr += 1;                   // Go to index result[3]
    *ptr = v[3] * s;


    return result;
}

/*  @function vecScalarMult_inplace:
 *      Calculate the Scalar multiplacation of the vector given, 
 *      and store the result in the third parameter.
 */
void vecScalarMult_inplace(GLfloat s, vec4 v, vec4 result) {
    result[0] = v[0] * s;
    result[1] = v[1] * s;
    result[2] = v[2] * s;
    result[3] = v[3] * s;
}

/*  @function vecAdd:
 *      Sum the two given vectors and return the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 */
vec4* vecAdd(vec4 vectorA, vec4 vectorB) {
    vec4* result = malloc(sizeof(vec4));
    GLfloat* ptr = result[0];

    *ptr = vectorA[0] + vectorB[0];
    ptr += 1;       // Go to index result[1]
    *ptr = vectorA[1] + vectorB[1];
    ptr += 1;       // Go to index result[2]
    *ptr = vectorA[2] + vectorB[2];
    ptr += 1;       // Go to index result[3]
    *ptr = vectorA[3] + vectorB[3];

    return result;
}

/*  @function vecScalarMult_inplace:
 *      Calculate the Scalar multiplacation of the vector given, 
 *      and store the result in the third parameter.
 */
void vecAdd_inplace(vec4 vectorA, vec4 vectorB, vec4 result) {
    result[0] = vectorA[0] + vectorB[0];
    result[1] = vectorA[1] + vectorB[1];
    result[2] = vectorA[2] + vectorB[2];
    result[3] = vectorA[3] + vectorB[3];
}

/*  @function vecSub:
 *      Subtract the two given vectors and return the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 */
vec4* vecSub(vec4 vectorA, vec4 vectorB) {
    vec4* result = malloc(sizeof(vec4));
    GLfloat* ptr = result[0];

    *ptr = vectorA[0] - vectorB[0];
    ptr += 1;
    *ptr = vectorA[1] - vectorB[1];
    ptr += 1;
    *ptr = vectorA[2] - vectorB[2];
    ptr += 1;
    *ptr = vectorA[3] - vectorB[3];

    return result;
}

/*  @function vecSub_inplace:
 *      Subtract the two given vectors and return the result in the third parameter (result)
 */
void vecSub_inplace(vec4 vectorA, vec4 vectorB, vec4 result) {
    result[0] = vectorA[0] - vectorB[0];
    result[1] = vectorA[1] - vectorB[1];
    result[2] = vectorA[2] - vectorB[2];
    result[3] = vectorA[3] - vectorB[3];
}

/*  @function vecDotProd:
 *      Calculate the vector dot product and return the resulting scalar. 
 */
GLfloat vecDotProd(vec4 vectorA, vec4 vectorB) {
    GLfloat result;
    result = vectorA[0] * vectorB[0];
    result += vectorA[1] * vectorB[1];
    result += vectorA[2] * vectorB[2];
    result += vectorA[3] * vectorB[3];

    return result;
}

/*  @function vecCrossProd:
 *      Calculate the vector cross product the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 */
vec4* vecCrossProd(vec4 vectorA, vec4 vectorB) {
    vec4* result = malloc(sizeof(vec4)); 
    GLfloat* ptr = result[0];   // Start at index result[0]

    *ptr = (vectorA[1] * vectorB[2]) - (vectorA[2] * vectorB[1]);   // result[0] = (y1 * z2) - (z1 * y2)
    ptr += 1;   // Go to index result[1]
    *ptr = (vectorA[2] * vectorB[0]) - (vectorA[0] * vectorB[2]);   // result[1] = (z1 * x2) - (x1 * z2)
    ptr += 1;   // Go to index result[2]
    *ptr = (vectorA[0] * vectorB[1]) - (vectorA[1] * vectorB[0]);   // result[2] = (x1 * y2) - (y1 * x2)
    ptr += 1;   // Go to index result[3]
    *ptr = 0;   // 4th element (w) should always be zero.           // result[3] = 0

    return result;
}

/*  @function vecCrossProd_inplace:
 *      Calculate the vector cross product the result in the 3rd parameter (result)
 */
void vecCrossProd_inplace(vec4 vectorA, vec4 vectorB, vec4 result) {
    result[0] = (vectorA[1] * vectorB[2]) - (vectorA[2] * vectorB[1]);  // result[0] = (y1 * z2) - (z1 * y2)
    result[1] = (vectorA[2] * vectorB[0]) - (vectorA[0] * vectorB[2]);  // result[1] = (z1 * x2) - (x1 * z2)
    result[2] = (vectorA[0] * vectorB[1]) - (vectorA[1] * vectorB[0]);  // result[2] = (x1 * y2) - (y1 * x2)
    result[3] = 0;                                                      // result[3] = 0
}
/**
 * @brief Get magnitude of a given 4x1 vector
 * 
 * @param v - The vector we want the magnitude of
 * @return GLfloat - The magnitude of the given vector
 */
GLfloat vecMagnitude(vec4 v) {
    GLfloat mag_2 = (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]);
    return sqrtf(mag_2);
}

/**
 * @brief Normalize the given 4x1 vector
 * 
 * @param v - The vector to normalize
 * @return vec4* - The normalized 4x1 vector
 */
vec4* vecNormalize(vec4 v) {
    GLfloat mag = vecMagnitude(v);
    mag = 1 / mag;
    return vecScalarMult(mag, v);;
}

/**
 * @brief Normalize the given 4x1 vector inplace
 * 
 * @param v - the vector to normalize
 * @param norm - the vector to store the normalized vector in
 */
void vecNormalize_inplace(vec4 v, vec4 norm) {
    GLfloat mag = vecMagnitude(v);
    mag = 1 / mag;
    vecScalarMult_inplace(mag, v, norm);
}

/**
 * @brief Get the vector object based off of two points
 *      Calculate the difference of two points
 * 
 * @param pointA - The first point
 * @param pointB - The seconnd point 
 * @return vec4* - The vector resulting from the two points
 */
vec4* get_vector(vec4 pointA, vec4 pointB) {
    return vecSub(pointB, pointA);
}

void get_vector_inplace(vec4 pointA, vec4 pointB, vec4 result) {
    vecSub_inplace(pointB, pointA, result);
}


/**
 * @brief Get the unit vector object based off of two points
 *      Calculate the difference of two points, and then normalize
 * 
 * @param pointA - The first point
 * @param pointB - The seconnd point 
 * @return vec4* - The unit vector resulting from the two points
 */
vec4* get_unit_vector(vec4 pointA, vec4 pointB) {
    vec4* u = get_vector(pointA, pointB);
    vec4* unit = vecNormalize(*u);
    free(u);
    return unit;
}

/***      Matrix Operations      ***/

/*  @function matPrint:
 *      Print each element of the 4x4 matrix with 3 points of floating
 *      point granularity
 */
void matPrint(mat4 mat) {
    printf("[%8.3f, %8.3f, %8.3f, %8.3f]\n"
           "[%8.3f, %8.3f, %8.3f, %8.3f]\n"
           "[%8.3f, %8.3f, %8.3f, %8.3f]\n"
           "[%8.3f, %8.3f, %8.3f, %8.3f]\n", 
           mat.x[0], mat.y[0], mat.z[0], mat.w[0],
           mat.x[1], mat.y[1], mat.z[1], mat.w[1],
           mat.x[2], mat.y[2], mat.z[2], mat.w[2],
           mat.x[3], mat.y[3], mat.z[3], mat.w[3]);
}

void mat3Print(mat3 mat) {
    printf("[%8.3f, %8.3f, %8.3f]\n"
           "[%8.3f, %8.3f, %8.3f]\n"
           "[%8.3f, %8.3f, %8.3f]\n", 
           mat.x[0], mat.y[0], mat.z[0],
           mat.x[1], mat.y[1], mat.z[1],
           mat.x[2], mat.y[2], mat.z[2]);
}

/*  @function matPrintRaw:
 *      Print each element of the 4x4 matrix with full floating
 *      point granularity (prints complete floating point value).
 */
void matPrintAll(mat4 mat) {
    printf("[%10f, %10f, %10f, %10f]\n"
           "[%10f, %10f, %10f, %10f]\n"
           "[%10f, %10f, %10f, %10f]\n"
           "[%10f, %10f, %10f, %10f]\n", 
           mat.x[0], mat.y[0], mat.z[0], mat.w[0],
           mat.x[1], mat.y[1], mat.z[1], mat.w[1],
           mat.x[2], mat.y[2], mat.z[2], mat.w[2],
           mat.x[3], mat.y[3], mat.z[3], mat.w[3]);
}

/*  @function matScalarMult:
 *      Compute scalar multiple of 4x4 matrix. 
 *      Multiply each column order vector by the given scalar value.
 *      returns 4x4 matris
 */
mat4* matScalarMult(GLfloat s, mat4* m) {
    mat4 *result = malloc(sizeof(mat4));

    vecScalarMult_inplace(s, m->x, result->x);
    vecScalarMult_inplace(s, m->y, result->y);
    vecScalarMult_inplace(s, m->z, result->z);
    vecScalarMult_inplace(s, m->w, result->w);

    return result;
}

/* Inplace version of matScalarMult
 *      Does not return anything, instead stores result in the parameter 'result'
 */
void matScalarMult_inplace(GLfloat s, mat4* m, mat4* result) {
    vecScalarMult_inplace(s, m->x, result->x);
    vecScalarMult_inplace(s, m->y, result->y);
    vecScalarMult_inplace(s, m->z, result->z);
    vecScalarMult_inplace(s, m->w, result->w);
}

/*  @function matAdd:
 *      Sum the two given 4x4 matrices column by column, and store result in third parameter. 
 */
mat4* matAdd(mat4* matrixA, mat4* matrixB) {
    mat4 *result = malloc(sizeof(mat4));

    vecAdd_inplace(matrixA->x, matrixB->x, result->x);
    vecAdd_inplace(matrixA->y, matrixB->y, result->y);
    vecAdd_inplace(matrixA->z, matrixB->z, result->z);
    vecAdd_inplace(matrixA->w, matrixB->w, result->w);

    return result;
}

/* Inplace version of matAdd
 *      Does not return anything, instead stores result in the parameter 'result'
 */
void matAdd_inplace(mat4* matrixA, mat4* matrixB, mat4* result) {
    vecAdd_inplace(matrixA->x, matrixB->x, result->x);
    vecAdd_inplace(matrixA->y, matrixB->y, result->y);
    vecAdd_inplace(matrixA->z, matrixB->z, result->z);
    vecAdd_inplace(matrixA->w, matrixB->w, result->w);
}

/*  @function matSub:
 *      Subtract the two given 4x4 matrices column by column, and store result in third parameter. 
 */
mat4* matSub(mat4* matrixA, mat4* matrixB) {
    mat4 *result = malloc(sizeof(mat4));

    vecSub_inplace(matrixA->x, matrixB->x, result->x);
    vecSub_inplace(matrixA->y, matrixB->y, result->y);
    vecSub_inplace(matrixA->z, matrixB->z, result->z);
    vecSub_inplace(matrixA->w, matrixB->w, result->w);

    return result;
}

/* Inplace version of matSub
 *      Does not return anything, instead stores result in the parameter 'result'
 */
void matSub_inplace(mat4* matrixA, mat4* matrixB, mat4* result) {
    vecSub_inplace(matrixA->x, matrixB->x, result->x);
    vecSub_inplace(matrixA->y, matrixB->y, result->y);
    vecSub_inplace(matrixA->z, matrixB->z, result->z);
    vecSub_inplace(matrixA->w, matrixB->w, result->w);
}

/*  @function matMult:
 *      Multiply the two given 4x4 matrices column by column, and store result in third parameter. 
 */
mat4* matMult(mat4* matrixA, mat4* matrixB) {
    mat4 *result = malloc(sizeof(mat4));

    vec4 matrixA_row_one = {matrixA->x[0], matrixA->y[0], matrixA->z[0], matrixA->w[0]};
    result->x[0] = vecDotProd(matrixA_row_one, matrixB->x);
    //printf("X IS MF : %f\n", result->x[0]);
    result->y[0] = vecDotProd(matrixA_row_one, matrixB->y);
    result->z[0] = vecDotProd(matrixA_row_one, matrixB->z);
    result->w[0] = vecDotProd(matrixA_row_one, matrixB->w);

    vec4 matrixA_row_two = {matrixA->x[1], matrixA->y[1], matrixA->z[1], matrixA->w[1]};
    result->x[1] = vecDotProd(matrixA_row_two, matrixB->x);
    result->y[1] = vecDotProd(matrixA_row_two, matrixB->y);
    result->z[1] = vecDotProd(matrixA_row_two, matrixB->z);
    result->w[1] = vecDotProd(matrixA_row_two, matrixB->w);

    vec4 matrixA_row_three = {matrixA->x[2], matrixA->y[2], matrixA->z[2], matrixA->w[2]};
    result->x[2] = vecDotProd(matrixA_row_three, matrixB->x);
    result->y[2] = vecDotProd(matrixA_row_three, matrixB->y);
    result->z[2] = vecDotProd(matrixA_row_three, matrixB->z);
    result->w[2] = vecDotProd(matrixA_row_three, matrixB->w);

    vec4 matrixA_row_four = {matrixA->x[3], matrixA->y[3], matrixA->z[3], matrixA->w[3]};
    result->x[3] = vecDotProd(matrixA_row_four, matrixB->x);
    result->y[3] = vecDotProd(matrixA_row_four, matrixB->y);
    result->z[3] = vecDotProd(matrixA_row_four, matrixB->z);
    result->w[3] = vecDotProd(matrixA_row_four, matrixB->w);

    return result;
}

/* Inplace version of matMult
 *      Does not return anything, instead stores result in the parameter 'result'
 */
void matMult_inplace(mat4* matrixA, mat4* matrixB, mat4* result) {
    vec4 matrixA_row_one = {matrixA->x[0], matrixA->y[0], matrixA->z[0], matrixA->w[0]};
    result->x[0] = vecDotProd(matrixA_row_one, matrixB->x);
    result->y[0] = vecDotProd(matrixA_row_one, matrixB->y);
    result->z[0] = vecDotProd(matrixA_row_one, matrixB->z);
    result->w[0] = vecDotProd(matrixA_row_one, matrixB->w);

    vec4 matrixA_row_two = {matrixA->x[1], matrixA->y[1], matrixA->z[1], matrixA->w[1]};
    result->x[1] = vecDotProd(matrixA_row_two, matrixB->x);
    result->y[1] = vecDotProd(matrixA_row_two, matrixB->y);
    result->z[1] = vecDotProd(matrixA_row_two, matrixB->z);
    result->w[1] = vecDotProd(matrixA_row_two, matrixB->w);

    vec4 matrixA_row_three = {matrixA->x[2], matrixA->y[2], matrixA->z[2], matrixA->w[2]};
    result->x[2] = vecDotProd(matrixA_row_three, matrixB->x);
    result->y[2] = vecDotProd(matrixA_row_three, matrixB->y);
    result->z[2] = vecDotProd(matrixA_row_three, matrixB->z);
    result->w[2] = vecDotProd(matrixA_row_three, matrixB->w);

    vec4 matrixA_row_four = {matrixA->x[3], matrixA->y[3], matrixA->z[3], matrixA->w[3]};
    result->x[3] = vecDotProd(matrixA_row_four, matrixB->x);
    result->y[3] = vecDotProd(matrixA_row_four, matrixB->y);
    result->z[3] = vecDotProd(matrixA_row_four, matrixB->z);
    result->w[3] = vecDotProd(matrixA_row_four, matrixB->w);
}

/*  @function matTranspose:
 *      Tranpose the give matrix and store into the 2nd parameter 
 */
mat4* matTranspose(mat4* m) {
    mat4 *result = malloc(sizeof(mat4));

    result->x[0] = m->x[0];
    result->x[1] = m->y[0];
    result->x[2] = m->z[0];
    result->x[3] = m->w[0];

    result->y[0] = m->x[1];
    result->y[1] = m->y[1];
    result->y[2] = m->z[1];
    result->y[3] = m->w[1];

    result->z[0] = m->x[2];
    result->z[1] = m->y[2];
    result->z[2] = m->z[2];
    result->z[3] = m->w[2];

    result->w[0] = m->x[3];
    result->w[1] = m->y[3];
    result->w[2] = m->z[3];
    result->w[3] = m->w[3];

    return result;
}

/* Inplace version of matTranspose
 *      Does not return anything, instead stores result in the parameter 'transpose'
 */
void matTranspose_inplace(mat4* m, mat4* transpose) {
    transpose->x[0] = m->x[0];
    transpose->x[1] = m->y[0];
    transpose->x[2] = m->z[0];
    transpose->x[3] = m->w[0];

    transpose->y[0] = m->x[1];
    transpose->y[1] = m->y[1];
    transpose->y[2] = m->z[1];
    transpose->y[3] = m->w[1];

    transpose->z[0] = m->x[2];
    transpose->z[1] = m->y[2];
    transpose->z[2] = m->z[2];
    transpose->z[3] = m->w[2];

    transpose->w[0] = m->x[3];
    transpose->w[1] = m->y[3];
    transpose->w[2] = m->z[3];
    transpose->w[3] = m->w[3];
}

/*  @function matVectorMult:
 *      Vector multiply matrix m by vec and store the resulting vector 
 */
vec4* matVectorMult(mat4* m, vec4 vec) {
    vec4 *result = malloc(sizeof(vec4));

    vec4 xColumn;
    vec4 yColumn;
    vec4 zColumn;
    vec4 wColumn;

    vecScalarMult_inplace(vec[0], m->x, xColumn);   // Multiply v[0] by matrix column x
    vecScalarMult_inplace(vec[1], m->y, yColumn);   // Multiply v[1] by matrix column y
    vecScalarMult_inplace(vec[2], m->z, zColumn);   // Multiply v[2] by matrix column z
    vecScalarMult_inplace(vec[3], m->w, wColumn);   // Multiply v[3] by matrix column w

    // Sum each row of 4x4 matrix and store into 4x1 vector
    GLfloat *ptr = result[0];   // Start at result[0]
    *ptr = xColumn[0] + yColumn[0] + zColumn[0] + wColumn[0];  
    ptr += 1;                   // Go to result[1]
    *ptr = xColumn[1] + yColumn[1] + zColumn[1] + wColumn[1];
    ptr += 1;                   // Go to result[2]
    *ptr = xColumn[2] + yColumn[2] + zColumn[2] + wColumn[2];
    ptr += 1;                   // Go to result[3]
    *ptr = xColumn[3] + yColumn[3] + zColumn[3] + wColumn[3];

    return result;
}

/* Inplace version of matVectorMult
 *      Does not return anything, instead stores result in the parameter 'result'
 */
void matVectorMult_inplace(mat4* m, vec4 vec, vec4 result) {
    vec4 xColumn;
    vec4 yColumn;
    vec4 zColumn;
    vec4 wColumn;

    vecScalarMult_inplace(vec[0], m->x, xColumn);   // Multiply v[0] by matrix column x
    vecScalarMult_inplace(vec[1], m->y, yColumn);   // Multiply v[1] by matrix column y
    vecScalarMult_inplace(vec[2], m->z, zColumn);   // Multiply v[2] by matrix column z
    vecScalarMult_inplace(vec[3], m->w, wColumn);   // Multiply v[3] by matrix column w

    // Sum each row of 4x4 matrix and store into 4x1 vector
    result[0] = xColumn[0] + yColumn[0] + zColumn[0] + wColumn[0];  
    result[1] = xColumn[1] + yColumn[1] + zColumn[1] + wColumn[1];
    result[2] = xColumn[2] + yColumn[2] + zColumn[2] + wColumn[2];
    result[3] = xColumn[3] + yColumn[3] + zColumn[3] + wColumn[3];
}

/********* Inverse (and all of it's dependent functions) ***********/ 

/*  @function matVectorMult:
 *      Calculates the inverse of matrix m. 
 *      If the inverse does not exist, it should return a 4x4 matrix of NaN
 */
mat4* matInverse(mat4* m) {
    mat4 *inverse = malloc(sizeof(mat4));

    mat4 minor;
    mat4 cofactor;
    mat4 transpose;
    GLfloat determinant;

    // Get the minor of matrix m, store in minor
    matMinors(m, &minor);

    // Get the cofactor of the matrix of minor, store in cofactor
    matCofactor(&minor, &cofactor);

    // Tranpose the matrix cofactor, store in transpose
    matTranspose_inplace(&cofactor, &transpose);

    // Get the determinant of m
    determinant = matDeterminant(m);

    // Final step to calculating the inverse
    matScalarMult_inplace(1/determinant, &transpose, inverse);

    return inverse;
}

/* Inplace version of matInverse
 *      Does not return anything, instead stores result in the parameter 'inverse'
 */
void matInverse_inplace(mat4* m, mat4* inverse) {
    mat4 minor;
    mat4 cofactor;
    mat4 transpose;
    GLfloat determinant;

    // Get the minor of matrix m, store in minor
    matMinors(m, &minor);

    // Get the cofactor of the matrix of minor, store in cofactor
    matCofactor(&minor, &cofactor);

    // Tranpose the matrix cofactor, store in transpose
    matTranspose_inplace(&cofactor, &transpose);

    // Get the determinant of m
    determinant = matDeterminant(m);

    // Final step to calculating the inverse
    matScalarMult_inplace(1/determinant, &transpose, inverse);
}

/*  @function matMinor:
 *      Get the minor of matrix m
 */  
void matMinors(mat4* m, mat4* minor) {
    mat3 subMatrix;

    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            // For each current element, reduce matrix to 3x3 by removing element's row and column
            mat4ToMat3(m, row, col, &subMatrix);
            if (col == 0) {         // Set X column 
                minor->x[row] = mat3Determinant(&subMatrix);
            } else if (col == 1) {  // Set Y column 
                minor->y[row] = mat3Determinant(&subMatrix);
            } else if (col == 2) {  // Set Z column 
                minor->z[row] = mat3Determinant(&subMatrix);
            } else if (col == 3) {  // Set W column 
                minor->w[row] = mat3Determinant(&subMatrix);
            }
        }
    }
}

/*  @function matCofactor:
 *      Get the cofactor of matrix m
 */  
void matCofactor(mat4* m, mat4* cofactor) {
    int a = 1;
    int b = -1;
    // Loop through each row making a checkerboard of [1 / (-1)]
    for (int row = 0; row < 4; row++) {
        cofactor->x[row] = m->x[row] * a;   // Multiply by a to get either positive or negative
        cofactor->y[row] = m->y[row] * b;
        cofactor->z[row] = m->z[row] * a;
        cofactor->w[row] = m->w[row] * b;

        a *= -1;    // Swap a +/-
        b *= -1;    // Swap b -/+
    }
}

/*  @function matDeterminant:
 *      Get the determinant of a 4x4 matrix.
 *      
 *      Break the 4x4 matrix into 4 3x3 matrices and get their determinants
 */ 
GLfloat matDeterminant(mat4* m) {
    mat3 m3;

    mat4ToMat3(m, 0, 0, &m3);   // Reduce 4x4 matrix into 3x3 matrix by removing xColumn and row 0
    GLfloat x = (m->x[0] * mat3Determinant(&m3));   // Calculate 3x3 determinant and multiply by xColumn Header

    mat4ToMat3(m, 0, 1, &m3);   // Reduce 4x4 matrix into 3x3 matrix by removing yColumn and row 0
    GLfloat y = (m->y[0] * mat3Determinant(&m3));   // Calculate 3x3 determinant and multiply by yColumn Header

    mat4ToMat3(m, 0, 2, &m3);   // Reduce 4x4 matrix into 3x3 matrix by removing zColumn and row 0
    GLfloat z = (m->z[0] * mat3Determinant(&m3));   // Calculate 3x3 determinant and multiply by zColumn Header

    mat4ToMat3(m, 0, 3, &m3);   // Reduce 4x4 matrix into 3x3 matrix by removing wColumn and row 0
    GLfloat w = (m->w[0] * mat3Determinant(&m3));   // Calculate 3x3 determinant and multiply by wColumn Header

    return (x - y + z - w);
}

/*  @function mat3Determinant:
 *      Return the determinant value of a 3x3 matrix
 */ 
GLfloat mat3Determinant(mat3* m) {
    GLfloat a = m->x[0] * ((m->y[1] * m->z[2]) - (m->z[1] * m->y[2]));
    GLfloat b = m->y[0] * ((m->x[1] * m->z[2]) - (m->x[2] * m->z[1]));
    GLfloat c = m->z[0] * ((m->x[1] * m->y[2]) - (m->x[2] * m->y[1]));

    return (a - b + c);
}

/**  @function mat4ToMat3:
 *      Reduce a 4x4 matrix into a 3x3 matrix
 *      Remove row i and column j and store the resulting 3x3 matrix in result.
 */ 
void mat4ToMat3(mat4* m4, int row, int col, mat3* result) {
    int currentRow = 0;
    for (int i = 0; i < 4; i++) {
        if (i != row) {
            if (col == 0) {         // Remove X column 
                result->x[currentRow] = m4->y[i];
                result->y[currentRow] = m4->z[i];
                result->z[currentRow] = m4->w[i];
                currentRow += 1;
            } else if (col == 1) {  // Remove Y column 
                result->x[currentRow] = m4->x[i];
                result->y[currentRow] = m4->z[i];
                result->z[currentRow] = m4->w[i];
                currentRow += 1;
            } else if (col == 2) {  // Remove Z column 
                result->x[currentRow] = m4->x[i];
                result->y[currentRow] = m4->y[i];
                result->z[currentRow] = m4->w[i];
                currentRow += 1;
            } else if (col == 3) {  // Remove W column 
                result->x[currentRow] = m4->x[i];
                result->y[currentRow] = m4->y[i];
                result->z[currentRow] = m4->z[i];
                currentRow += 1;
            }
        }
    }
}


/*****************************************************************************
 * 
 *              Extended Functionality: Rotate/Translate/Scale/etc
 *
 *****************************************************************************/
/**
 * @brief Concatenate/"chain" together a variable number of tranformations.
 * 
 * @param num - The number of tranformations you want to concatenate ('chain') together
 * @param ... - Variable argument for the tranformations you want to apply to the 
 *              Expects a pointer of mat4
 * @return mat4* returns the tranformation matrix (ctm) after concatenating all the given tranformations
 */
mat4* chainCTM(int num, ...) {
    // Allocate space for affine transformation matrix
    mat4* ctm = (mat4 *)calloc(1, sizeof(mat4));

    // Set ctm as identity matrix
    ctm->x[0] = 1;
    ctm->y[1] = 1;
    ctm->z[2] = 1;
    ctm->w[3] = 1;

    va_list args;           // List of arguments
    va_start(args, num);    // Initialize args as as an argument list with size = num
    
    mat4* currMatrix;

    for (int i = 0; i < num; i++) {
        currMatrix = va_arg(args, mat4*);

        ctm = matMult(ctm, currMatrix);

        free(currMatrix);
    }
    va_end(args);

    return ctm;
}

/**
 * @brief Concatenate/"chain" together a variable number of tranformations.
 * 
 * @param num - The number of tranformations you want to concatenate ('chain') together
 * @param ... - Variable argument for the tranformations you want to apply to the 
 *              Expects a pointer of mat4
 * @return mat4* returns the tranformation matrix (ctm) after concatenating all the given tranformations
 */
// TODO Fix the order of transformations
mat4* chainCurrCTM(mat4 *ctm, int num, ...) {
    va_list args;           // List of arguments
    va_start(args, num);    // Initialize args as as an argument list with size = num

    mat4* currMatrix;

    for (int i = 0; i < num; i++) {
        currMatrix = va_arg(args, mat4*);

        ctm = matMult(ctm, currMatrix);

        free(currMatrix);
    }
    va_end(args);

    return ctm;
}

/**
 * @brief Get a tranformation matrix that will translate the vertices by
 *          x, y, z.
 * 
 * @param x - translates x-axis
 * @param y - translate on y-axis by y
 * @param z - translate on z-axis by z
 * @return mat4* - The tranformation matrix of the translation.
 */
mat4* translate(GLfloat x, GLfloat y, GLfloat z) {
    // Remember, the mat4 is in column major order.
    mat4 *translationMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Basically make an indentity matrix    
    translationMatrix->x[0] = 1;
    translationMatrix->y[1] = 1;
    translationMatrix->z[2] = 1;
    translationMatrix->w[3] = 1;

    // Set the values x, y, and z values for the tranlation matrix
    translationMatrix->w[0] = x;
    translationMatrix->w[1] = y;
    translationMatrix->w[2] = z;
    
    // Return pointPrime (vec4)
    return translationMatrix;
}

mat4* translate_v(vec4 v) {
    return translate(v[0], v[1], v[2]);
}

/**
 * @brief Create a scaling matrix by the given scaling values. 
 * 
 * @param xScale - Value to scale the x-axis.
 * @param yScale - Value to scale the y-axis
 * @param zScale - Value to scale the z-axis
 * @return mat4* - Returns transformation matrix of the scaling
 */
mat4* scale(GLfloat xScale, GLfloat yScale, GLfloat zScale) {
    mat4 *scaleMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Set the scales for x, y, and z
    scaleMatrix->x[0] = xScale;
    scaleMatrix->y[1] = yScale;
    scaleMatrix->z[2] = zScale;
    scaleMatrix->w[3] = 1;

    // Return pointPrime (vec4)
    return scaleMatrix;
}

/**
 * @brief Creates an inverted scaling matrix. 
 *        Essentially reverts the scaling.
 *        
 *        scale_inv(a, b, c) * scale(a, b, c) => No change; i.e they cancel each other out.
 *        
 * 
 * @param xScale - Value to scale the x-axis; Sets as 1 / xScale
 * @param yScale - Value to scale the y-axis; Sets as 1 / yScale
 * @param zScale - Value to scale the z-axis; Sets as 1 / zScale
 * @return mat4* - Returns inverted scaling tranformation
 */
mat4* scale_inv(GLfloat xScale, GLfloat yScale, GLfloat zScale) {
    mat4 *scaleMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Set the scales for x, y, and z
    scaleMatrix->x[0] = 1.0f / xScale;
    scaleMatrix->y[1] = 1.0f / yScale;
    scaleMatrix->z[2] = 1.0f / zScale;
    scaleMatrix->w[3] = 1;

    // Return point of pointPrime (vec4)
    return scaleMatrix;
}

/**
 * @brief Rotate about the x-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param deg - The degree angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_x(GLfloat deg) {
    GLfloat rad = (deg * M_PI) / 180;

    return rotate_x_rad(rad);
}

/**
 * @brief Rotate about the x-axis by radians
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_x_rad(GLfloat rad) {
    mat4 *rotationMatrix = (mat4 *)calloc(1, sizeof(mat4));
    
    rotationMatrix->x[0] = 1;

    // y cos(rad) - z sin(rad)
    rotationMatrix->y[1] = cosf(rad);
    rotationMatrix->z[1] = -sinf(rad);

    // y sin(rad) + z cos(rad)
    rotationMatrix->y[2] = sinf(rad);
    rotationMatrix->z[2] = cosf(rad);

    rotationMatrix->w[3] = 1;

    return rotationMatrix;
}

/**
 * @brief Rotate to the x-axis using the parameters 
 *
 * @param alpha_y - The y value of the axis's of rotation normalized vector value
 * @param alpha_z - The z value of the axis's of rotation normalized vector value.
 * @param d - the value from the sqrt(y^2 + z^2)
 *              y and z are values of the normalized vector axis of rotation
 * @return mat4* - Returns the rotation matrix that will rotate the axis vector to the y-axis.
 */
mat4* rotate_to_x(GLfloat alpha_y, GLfloat alpha_z, GLfloat d) {
    mat4 *rotationMatrix = (mat4 *)calloc(1, sizeof(mat4));
    
    /* Create Identity Matrix */
    rotationMatrix->x[0] = 1;
    rotationMatrix->y[1] = 1;
    rotationMatrix->z[2] = 1;
    rotationMatrix->w[3] = 1;

    if (d != 0) {
        // 0 | a_z/d | -a_y/d | 0
        rotationMatrix->y[1] = alpha_z/d;
        rotationMatrix->z[1] = -alpha_y/d;

        // 0 | a_y/d | a_z/d | 0
        rotationMatrix->y[2] = alpha_y/d;
        rotationMatrix->z[2] = alpha_z/d;
    }

    return rotationMatrix;
}

/**
 * @brief Rotate about the y-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param deg - The degree angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_y(GLfloat deg) {
    GLfloat rad = (deg * M_PI) / 180;

    return rotate_y_rad(rad);
}

/**
 * @brief Rotate about the y-axis by degrees.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_y_rad(GLfloat rad) {
    mat4 *rotationMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // z sin(rad) + x cos(rad)
    rotationMatrix->z[0] = sinf(rad);
    rotationMatrix->x[0] = cosf(rad);

    rotationMatrix->y[1] = 1;

    // z cos(rad) - x sin(rad)
    rotationMatrix->z[2] = cosf(rad);
    rotationMatrix->x[2] = -sinf(rad);
    
    rotationMatrix->w[3] = 1;

    return rotationMatrix;
}

/**
 * @brief Rotate to the y-axis using the parameters 
 * 
 * @param alpha_x - The x value of the axis's of rotation normalized vector value
 * @param d - the value from the sqrt(y^2 + z^2)
 *              y and z are values of the normalized vector axis of rotation
 * @return mat4* - Returns the rotation matrix that will rotate the axis vector to the y-axis.
 */
mat4* rotate_to_y(GLfloat alpha_x, GLfloat d) {
    mat4 *rotationMatrix = (mat4 *)calloc(1, sizeof(mat4));

    /* Create Identity Matrix */
    rotationMatrix->x[0] = 1;
    rotationMatrix->y[1] = 1;
    rotationMatrix->z[2] = 1;
    rotationMatrix->w[3] = 1;

    // z sin(rad) + x cos(rad)
    if (d != 0) {
        rotationMatrix->x[0] = d;
        rotationMatrix->z[0] = -alpha_x;

        //rotationMatrix->y[1] = 1;

        // z cos(rad) - x sin(rad)
        rotationMatrix->x[2] = alpha_x;
        rotationMatrix->z[2] = d;
        
        //rotationMatrix->w[3] = 1;
    } else {
        rotationMatrix->x[0] = 0;
        rotationMatrix->x[2] = alpha_x;
        rotationMatrix->z[0] = -alpha_x;
        rotationMatrix->z[2] = 0;
    }

    return rotationMatrix;
}

/**
 * @brief Rotate about the z-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_z(GLfloat deg) {
    GLfloat rad = (deg * M_PI) / 180;

    return rotate_z_rad(rad);
}

/**
 * @brief Rotate about the z-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_z_rad(GLfloat rad) {
    mat4 *rotationMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // x cos(rad) - y sin(rad)
    rotationMatrix->x[0] = cosf(rad);
    rotationMatrix->y[0] = -sinf(rad);

    // x sin(rad) + y cos(rad)
    rotationMatrix->x[1] = sinf(rad);
    rotationMatrix->y[1] = cosf(rad);

    rotationMatrix->z[2] = 1;
    rotationMatrix->w[3] = 1;

    

    return rotationMatrix;
}

/**
 * @brief Shears about the x-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in degrees to shear about the y-axis
 * @param phi - The angle in degrees to shear about the z-axis
 * @return mat4* - The shearing matrix about the x-axis
 */
mat4* shear_x(GLfloat theta, GLfloat phi) {
    // Convert to radians
    theta = (theta * M_PI) / 180;
    phi = (phi * M_PI) / 180;

    return shear_x_rad(theta, phi);
}

/**
 * @brief Shears about the x-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in radians to shear about the y-axis
 * @param phi - The angle in radians to shear about the z-axis
 * @return mat4* - The shearing matrix about the x-axis
 */
mat4* shear_x_rad(GLfloat theta, GLfloat phi) {
    mat4 *shearMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Check theta and pi arent a multiple of 180. Trying to avoid dividing by zero error
    // tan(180) = 0
    // cot = 1/tan(180) -> Divide by zero error.
    if (fmod(theta, 180) == 0) { 
        theta = M_PI_4; // Set theta to pi/4 -> 90 degrees
    }
    if (fmod(phi, 180) == 0) {
        phi = M_PI_4; // Set phi to pi/4
    }

    shearMatrix->x[0] = 1;
    shearMatrix->y[1] = 1;
    shearMatrix->z[2] = 1;
    shearMatrix->w[3] = 1;

    shearMatrix->x[1] = 1 / tanf(theta);
    shearMatrix->x[2] = 1 / tanf(phi);

    return shearMatrix;
}

/**
 * @brief Shears about the y-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in degrees to shear about the x-axis
 * @param phi - The angle in degrees to shear about the z-axis
 * @return mat4* - The shearing matrix about the y-axis
 */
mat4* shear_y(GLfloat theta, GLfloat phi) {
    // Convert to radians
    theta = (theta * M_PI) / 180;
    phi = (phi * M_PI) / 180;

    return shear_y_rad(theta, phi);
}

/**
 * @brief Shears about the y-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in radians to shear about the x-axis
 * @param phi - The angle in radians to shear about the z-axis
 * @return mat4* - The shearing matrix about the y-axis
 */
mat4* shear_y_rad(GLfloat theta, GLfloat phi) {
    mat4 *shearMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Check theta and pi arent a multiple of 180. Trying to avoid dividing by zero error
    // tan(180) = 0
    // cot = 1/tan(180) -> Divide by zero error.
    if (fmod(theta, 180) == 0) { 
        theta = M_PI_4; // Set theta to pi/4 -> 90 degrees
    }
    if (fmod(phi, 180) == 0) {
        phi = M_PI_4; // Set phi to pi/4
    }

    shearMatrix->x[0] = 1;
    shearMatrix->y[1] = 1;
    shearMatrix->z[2] = 1;
    shearMatrix->w[3] = 1;
    
    shearMatrix->y[0] = 1 / tanf(theta);
    shearMatrix->y[2] = 1 / tanf(phi);

    return shearMatrix;
}

/**
 * @brief Shears about the z-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in degrees to shear about the x-axis
 * @param phi - The angle in degrees to shear about the y-axis
 * @return mat4* - The shearing matrix about the z-axis
 */
mat4* shear_z(GLfloat theta, GLfloat phi) {
    // Convert to radians
    theta = (theta * M_PI) / 180;
    phi = (phi * M_PI) / 180;

    return shear_z_rad(theta, phi);
}

/**
 * @brief Shears about the z-axis with and the center of shearing at the origin
 * 
 * @param theta - The angle in radians to shear about the x-axis
 * @param phi - The angle in radians to shear about the y-axis
 * @return mat4* - The shearing matrix about the z-axis
 */
mat4* shear_z_rad(GLfloat theta, GLfloat phi) {
    mat4 *shearMatrix = (mat4 *)calloc(1, sizeof(mat4));

    // Check theta and pi arent a multiple of 180. Trying to avoid dividing by zero error
    // tan(180) = 0
    // cot = 1/tan(180) -> Divide by zero error.
    if (fmod(theta, 180) == 0) { 
        theta = M_PI_4; // Set theta to pi/4 -> 90 degrees
    }
    if (fmod(phi, 180) == 0) {
        phi = M_PI_4; // Set phi to pi/4
    }

    shearMatrix->x[0] = 1;
    shearMatrix->y[1] = 1;
    shearMatrix->z[2] = 1;
    shearMatrix->w[3] = 1;

    shearMatrix->z[0] = 1 / tanf(theta);
    shearMatrix->z[1] = 1 / tanf(phi);

    return shearMatrix;
}

/**
 * @brief Rotate about a fixed point (point) about the x-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in degrees to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_x(vec4 point, GLfloat theta) {
    //Convert theta to radians
    theta = (theta * M_PI) / 180;

    return rotate_point_x_rad(point, theta);
}

/**
 * @brief Rotate about a fixed point (point) about the x-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in radians to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_x_rad(vec4 point, GLfloat theta) {
    return chainCTM(3,  translate(point[0], point[1], point[2]),
                        rotate_x_rad(theta),
                        translate(-point[0], -point[1], -point[2]));
}

/**
 * @brief Rotate about a fixed point (point) about the y-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in degrees to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_y(vec4 point, GLfloat theta) {
    // Convert theta to radians
    theta = (theta * M_PI) / 180;

    return rotate_point_y_rad(point, theta); 
}

/**
 * @brief Rotate about a fixed point (point) about the y-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in radians to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_y_rad(vec4 point, GLfloat theta) {
    return chainCTM(3,  translate(point[0], point[1], point[2]),
                        rotate_y_rad(theta),
                        translate(-point[0], -point[1], -point[2]));
}

/**
 * @brief Rotate about a fixed point (point) about the z-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in degrees to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_z(vec4 point, GLfloat theta) {
    theta = (theta * M_PI) / 180;

    return rotate_point_z_rad(point, theta);
}

/**
 * @brief Rotate about a fixed point (point) about the z-axis by the angle of ratation (theta)
 * 
 * @param point - The fixed point to ratote about
 * @param theta - The angle in radians to rotate by
 * @return mat4* - The rotation matrix about a fixed point 
 */
mat4* rotate_point_z_rad(vec4 point, GLfloat theta) {
    // Variadic arguments go from right to left.
    return chainCTM(3,  translate(point[0], point[1], point[2]),
                        rotate_z_rad(theta),
                        translate(-point[0], -point[1], -point[2]));
}

/**
 * @brief Rotate about an arbitrary axis by the angle of rotation (theta_z)
 * 
 * @param axis - The arbitrary axis to rotate about
 * @param theta_z - The angle in degrees to rotate by
 * @return mat4* - The rotation matrix about an arbitrary axis.
 */
mat4* rotate_axis(vec4 axis, GLfloat theta_z) {
    theta_z = (theta_z * M_PI) / 180;

    return rotate_axis_rad(axis, theta_z);
}

/**
 * @brief Rotate about an arbitrary axis by the angle of rotation (theta_z)
 * 
 * @param axis - The arbitrary axis to rotate about
 * @param theta_z - The angle in radians to rotate by
 * @return mat4* - The rotation matrix about an arbitrary axis.
 */
mat4* rotate_axis_rad(vec4 axis, GLfloat theta_z) { 
    mat4 *ctm;
    GLfloat a_x = axis[0];
    GLfloat a_y = axis[1];
    GLfloat a_z = axis[2];

    if (!isnan(a_x)) {
        GLfloat d = sqrtf((a_y * a_y) + (a_z * a_z));

        if (isnan(d)) 
            d = 0;

        // Calculate the transformationn matrices
        ctm = chainCTM(5, rotate_to_x(-a_y, a_z, d),
                                rotate_to_y(-a_x, d),
                                rotate_z_rad(theta_z),
                                rotate_to_y(a_x, d),
                                rotate_to_x(a_y, a_z, d));
    } else {
        ctm = calloc(1, sizeof(mat4));

        ctm->x[0] = 1;
        ctm->y[1] = 1;
        ctm->z[2] = 1;
        ctm->w[3] = 1;
    }

    return ctm;
}

/**
 * @brief Rotate about an arbitrary axis on a fixed point by the angle of rotation (theta_z)
 * 
 * @param fp - The fixed point to rotaate about
 * @param axis - The arbitrary axis to rotate about
 * @param theta_z - The angle in degrees to rotate by
 * @return mat4* - The rotation matrix about an arbitrary axis on a fixed point.
 */
mat4* rotate_point_axis(vec4 fp, vec4 axis, GLfloat theta_z) {
    theta_z = (theta_z * M_PI) / 180;

    return rotate_point_axis_rad(fp, axis, theta_z);
}

/**
 * @brief Rotate about an arbitrary axis on a fixed point by the angle of rotation (theta_z)
 * 
 * @param fp - The fixed point to rotate about
 * @param axis - The arbitrary axis to rotate about
 * @param theta_z - The angle in radians to rotate by
 * @return mat4* - The rotation matrix about an arbitrary axis on a fixed point.
 */
mat4* rotate_point_axis_rad(vec4 fp, vec4 axis, GLfloat theta_z) {
    return chainCTM(3,  translate(fp[0], fp[1], fp[2]),
                        rotate_axis_rad(axis, theta_z),
                        translate(-fp[0], -fp[1], -fp[2]));
}

