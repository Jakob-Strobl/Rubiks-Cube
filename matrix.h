

#ifndef MATRIX_H_
#define MATRIX_H_

#ifdef __APPLE__  // include Mac OS X verions of headers
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else // non-Mac OS X operating systems
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#endif  // __APPLE__

/***    Extern      ***/
#define identity_m4m4_matrix {{1, 0, 0, 0},{0, 1, 0, 0},{0, 0, 1, 0},{0, 0, 0, 1}}

/***    Typedef     ***/

/*  Define a type named 'vec4' that is a array of 4 Glfloats   
 *  @def vec4:     
 *      A simple 4x1 vector of GLfloat (typedef float) with,
 *          vec4[0] -> X axis   (Red)
 *          vec4[1] -> Y axis   (Green)
 *          vec4[2] -> Z axis   (Blue)
 *          vec4[3] -> W value  (Alpha)
 */
typedef GLfloat vec4[4];

/*  Define a type named 'vec3' that is a array of 3 Glfloats   
 *  @def vec3:     
 *      A simple 3x1 vector of GLfloat (typedef float) with,
 *          vec4[0] -> X axis 
 *          vec4[1] -> Y axis  
 *          vec4[2] -> Z axis  
 *  This is mainly used for getting the determinants and matrix of minor 
 */
typedef GLfloat vec3[3];

/*  Define a type named 'vec2' that is a array of 2 Glfloats   
 *  @def vec2:     
 *      A simple 3x1 vector of GLfloat (typedef float) with,
 *          vec4[0] -> X axis 
 *          vec4[1] -> Y axis   
 *  This is mainly used for textures and other types of 2D data
 */
typedef GLfloat vec2[2];

/*  Define a type named 'mat4' that is a struct of 4 vec4's (4x1 vectors)   
 *  @def mat4:     
 *      A 4x4 matrix composed of 4 column-major 4x1 vectors (type vec4)
 *          mat4.x -> X column 
 *          mat4.y -> Y column
 *          mat4.z -> Z column
 *          mat4.w -> W column 
 */
typedef struct {
    vec4 x;     /* X column vector (@see vec4 definition above) */
    vec4 y;     /* Y column vector  */
    vec4 z;     /* Z column vector */
    vec4 w;     /* w column vector */
} mat4;

/*  Define a type named 'mat3' that is a struct of 3 vec3's (3x1 vectors)   
 *  @def mat3:     
 *      A 3x3 matrix composed of 4 column-major 4x1 vectors (type vec4)
 *          mat4.x -> X column 
 *          mat4.y -> Y column
 *          mat4.z -> Z column
 * 
 *  This is mainly used for getting the determinants and matrix of minor 
 */
typedef struct {
    vec3 x;     /* X column vector (@see vec3 definition above) */
    vec3 y;     /* Y column vector  */
    vec3 z;     /* Z column vector  */
} mat3;

/*  Vector Operations */
/*  @function vecPrint:
 *      Print each element of the vector with 3 points of floating
 *      point granularity
 * 
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param vec4 vec - 4x1 vector of 4 GLfloats
 */
void vecPrint(vec4 v);

/*  @function vecPrintRaw:
 *      Print each element of the vector with full floating
 *      point granularity (prints complete floating point value).
 * 
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param vec4 vec - 4x1 vector of 4 GLfloats
 */
void vecPrintAll(vec4 v);

/*  @function vecScalarMult:
 *      Calculate the Scalar multiplacation of the vector given, and return malloced vec4
 *      Remember to free() the vector when finished. 
 * 
 *  @param vec4 vec - 4x1 vector of 4 GLfloats
 *  @param GLfloat s - scalar float value to multiply each element in the vector by. 
 *  @return vec4* - pointer to malloc'd vec4 (Glfloat[4])
 */
vec4* vecScalarMult(GLfloat s, vec4 v);

/*  @function vecScalarMult_inplace:
 *      Calculate the Scalar multiplacation of the vector given, 
 *      and store the result in the third parameter.
 * 
 *  @param vec4 vec - 4x1 vector of 4 GLfloats
 *  @param GLfloat s - scalar float value to multiply each element in the vector by. 
 *  @param result - stores result of vector scalar multiplication.
 */
void vecScalarMult_inplace(GLfloat s, vec4 v, vec4 result);

/*  @function vecAdd:
 *      Sum the two given vectors and return the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 * 
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @return vec4* - pointer to malloc'd vec4 (Glfloat[4])
 */
vec4* vecAdd(vec4 vectorA, vec4 vectorB);

/*  @function vecAdd_inplace:
 *      Sum the two given vectors and return the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @param vec4 result  - store vector addition
 */
void vecAdd_inplace(vec4 vectorA, vec4 vectorB, vec4 result);

/*  @function vecSub:
 *      Subtract the two given vectors and return the result in a new malloc'd vec4.
 *      Return = vA - vB 
 *      Remember to free() the vector when finished. 
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @return vec4* - pointer to malloc'd vec4 (Glfloat[4])
 */
vec4* vecSub(vec4 vectorA, vec4 vectorB);

/*  @function vecSub_inplace:
 *      Subtract the two given vectors and return the result in the third parameter (result)
 *      vA - vB = result
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @param vec4 result  - store vector addition
 */
void vecSub_inplace(vec4 vectorA, vec4 vectorB, vec4 result);

/*  @function vecDotProd:
 *      Calculate the vector dot product and return the resulting scalar. 
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @return GLfloat - the scalar result
 */
GLfloat vecDotProd(vec4 vectorA, vec4 vectorB);

/*  @function vecCrossProd:
 *      Calculate the vector cross product the result in a new malloc'd vec4.
 *      Remember to free() the vector when finished. 
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @return vec4* - the resulting cross product pointer to malloc'd vec4 (Glfloat[4])
 */
vec4* vecCrossProd(vec4 vectorA, vec4 vectorB);

/*  @function vecCrossProd_inplace:
 *      Calculate the vector cross product the result in the 3rd parameter (result)
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 *  @param vec4 vectorB - 4x1 vector of 4 GLfloats  
 *  @param vec4 result  - store vector cross product
 *  @return vec4* - the resulting cross product pointer to malloc'd vec4 (Glfloat[4])
 */
void vecCrossProd_inplace(vec4 vectorA, vec4 vectorB, vec4 result);

/*  @function vecMagnitude:
 *      Calculate the magnitude of a vector
 * 
 *  @param vec4 vectorA - 4x1 vector of 4 GLfloats
 */
GLfloat vecMagnitude(vec4 v);

/*  @function vecNormalize:
 *      Calculate the normalized vector
 * 
 *  @param vec4 vectorA - vector to be normalized
 *  @return vec4 - normalized vector.
 */
vec4* vecNormalize(vec4 v);

/*  @function vecNormalize_inplace:
 *      Calculate the normalized vector
 * 
 *  @param vec4 vectorA - vector to be normalized
 *  @param vec4 vectorA - normalized vector.
 */
void vecNormalize_inplace(vec4 v, vec4 norm);

vec4* get_vector(vec4 pointA, vec4 pointB);
void get_vector_inplace(vec4 pointA, vec4 pointB, vec4 result);
vec4* get_unit_vector(vec4 pointA, vec4 pointB);

/************************  
 *   Matrix Operations  *
 ************************/

/*  @function mat3Print:
 *      Print each element of the 3x3 matrix with 3 points of floating
 *      point granularity
 * 
 *  @see matrix.h for type-definitions of mat3 and vec3
 * 
 *  @param mat3 mat - 3x3 matrix composed of 3 column major 3x1 vectors (vec3)
 */
void mat3Print(mat3 m);


/*  @function matPrint:
 *      Print each element of the 4x4 matrix with 3 points of floating
 *      point granularity
 * 
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param mat4 mat - 4x4 matrix composed of 4 column major 4x1 vectors (vec4)
 */
void matPrint(mat4 m);

/*  @function matPrintAll:
 *      Print each element of the 4x4 matrix with full floating
 *      point granularity (prints complete floating point value).
 *  
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param mat4 mat - 4x4 matrix composed of 4 column major 4x1 vectors (vec4)
 */
void matPrintAll(mat4 m);

/*  @function matScalarMult:
 *      Compute scalar multiple of 4x4 matrix. 
 *      Multiply each column order vector by the given scalar value.
 *  
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  @param mat4 mat - 4x4 matrix composed of 4 column major 4x1 vectors (vec4)
 *  @param GLfloat s - scalar float value to multiply each element in the vector by. 
 *  @return mat4 - 4x4 matrix struct composed of 4 column order vec4s
 */
mat4* matScalarMult(GLfloat s, mat4* m);
void matScalarMult_inplace(GLfloat s, mat4* m, mat4* result);

/*  @function matAdd:
 *      Sum the two given 4x4 matrices column by column, and store result in third parameter. 
 *  
 *  @see matrix.h for type-definitions of mat4 and vec4
 * 
 *  mat4* - Pointer to 4x4 matrix of 4 column major 4x1 vectors (vec4)
 * 
 *  @param mat4* matrixA - Left-hand operator of summation 
 *  @param mat4* matrixB - Right-hand operator of summation
 *  @param mat4* result - Pointer of struct where the summation results will be stored (overwrites previous values)
 */
mat4* matAdd(mat4* matrixA, mat4* matrixB);
void matAdd_inplace(mat4* matrixA, mat4* matrixB, mat4* result);

/*  @function matSub:
 *      Subtract the two given 4x4 matrices column by column, and store result in third parameter. 
 *  
 *  @param mat4* matrixA - Left-hand operator of subtraction 
 *  @param mat4* matrixB - Right-hand operator of subtraction
 *  @param mat4* result - Pointer of struct where the subtraction results will be stored (overwrites previous values)
 */
mat4* matSub(mat4* matrixA, mat4* matrixB);
void matSub_inplace(mat4* matrixA, mat4* matrixB, mat4* result);

/*  @function matMult:
 *      Multiply the two given 4x4 matrices column by column, and store result in third parameter. 
 *  
 *  @param mat4* matrixA - Left-hand operator of multiplication 
 *  @param mat4* matrixB - Right-hand operator of multiplication
 *  @param mat4* result - Pointer of struct where the multiplication results will be stored (overwrites previous values)
 */
mat4* matMult(mat4* matrixA, mat4* matrixB);
void matMult_inplace(mat4* matrixA, mat4* matrixB, mat4* result);

/*  @function matTranspose:
 *      Transpose matrix m and store tranposed matrix (mT) in 2nd parameter (transpose)
 *  
 *  @param mat4* matrixA - Left-hand side of operation
 *  @param mat4* result - the resulting transposed matrix
 */
mat4* matTranspose(mat4* m);
void matTranspose_inplace(mat4* m, mat4* transpose);

/*  @function matVectorMult:
 *      Vector multiply matrix m by vec and store the resulting vector 
 *  
 *  @param mat4* m - Left-hand side of operation
 *  @param vec4* vec - Right-hand side of operation
 *  @param vec4 result -  result of vector multiplication
 */
vec4* matVectorMult(mat4* m, vec4 vec);
void matVectorMult_inplace(mat4* m, vec4 vec, vec4 result);

/*  @function matInverse:
 *      Calculate the inverse matrix of m
 *  
 *  @param mat4* m - Left-hand side of operation
 *  @param mat4* inverse -  the inverse of m
 */
mat4* matInverse(mat4* m);
void matInverse_inplace(mat4* m, mat4* inverse);
void matMinors(mat4* m, mat4* minor);
void matCofactor(mat4* m, mat4* cofactor);
GLfloat matDeterminant(mat4* m);
GLfloat mat3Determinant(mat3* m);

/**
 * @brief Convert a 4x4 matrix (mat4) into a 3x3 matrix (mat3)
 *          We choose a column and row to remove from a given 4x4 matrix.
 * 
 * @param m4 - The 4x4 matrix we want to convert into a 3x3 matrix
 * @param row - The row we want remove; -> get 4x3 matrix
 * @param col - The col we will remove; -> get 3x3 matrix
 * @param result - The 3x3 matrix.
 */
void mat4ToMat3(mat4* m4, int row, int col, mat3* result);


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
mat4* chainCTM(int num, ...);

mat4* chainCurrCTM(mat4 *ctm, int num, ...);

/**
 * @brief Get a tranformation matrix that will translate the vertices by
 *          x, y, z.
 * 
 * @param x - translates x-axis
 * @param y - translate on y-axis by y
 * @param z - translate on z-axis by z
 * @return mat4* - The tranformation matrix of the translation.
 */
mat4* translate(GLfloat x, GLfloat y, GLfloat z);

mat4* translate_v(vec4 v);

/**
 * @brief Create a scaling matrix by the given scaling values. 
 * 
 * @param xScale - Value to scale the x-axis.
 * @param yScale - Value to scale the y-axis
 * @param zScale - Value to scale the z-axis
 * @return mat4* - Returns transformation matrix of the scaling
 */
mat4* scale(GLfloat xScale, GLfloat yScale, GLfloat zScale);

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
mat4* scale_inv(GLfloat xScale, GLfloat yScale, GLfloat zScale);

/**
 * @brief Rotate around x-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param deg - The degree angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_x(GLfloat deg);

/**
 * @brief Rotate around x-axis by radians
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_x_rad(GLfloat rad);

/**
 * @brief Rotate around y-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param deg - The degree angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_y(GLfloat deg);

/**
 * @brief Rotate around y-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_y_rad(GLfloat rad);

/**
 * @brief Rotate around z-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_z(GLfloat deg);

/**
 * @brief Rotate around z-axis by degrees.
 *        Just converts degrees into radians, 
 *        and then invokes the radian function.
 * 
 * @param rad - The radian angle to rotate by
 * @return mat4* - Rotation tranformation matrix
 */
mat4* rotate_z_rad(GLfloat rad);

/**
 * @brief Shear along the x-axis
 * 
 * @param theta - Angle of shearing about y-axis.
 * @param phi - Angle of shearing about z-axis.
 *      Theta and Phi if = 90 degrees. Do nothing.
 * @return mat4* - Shearing matrix
 */
mat4* shear_x(GLfloat theta, GLfloat phi);
mat4* shear_x_rad(GLfloat theta, GLfloat phi);

/**
 * @brief Shear along the y-axis
 * 
 * @param theta - Angle of shearing about x-axis.
 * @param phi - Angle of shearing about z-axis.
 *      Theta and Phi if = 90 degrees. Do nothing.
 * @return mat4* - Shearing matrix
 */
mat4* shear_y(GLfloat theta, GLfloat phi);
mat4* shear_y_rad(GLfloat theta, GLfloat phi);

/**
 * @brief Shear along the z-axis
 * 
 * @param theta - Angle of shearing about x-axis.
 * @param phi - Angle of shearing about y-axis.
 *      Theta and Phi if = 90 degrees. Do nothing.
 * @return mat4* - Shearing matrix
 */
mat4* shear_z(GLfloat theta, GLfloat phi);
mat4* shear_z_rad(GLfloat theta, GLfloat phi);

/**
 * @brief 
 * 
 * @param point 
 * @param theta 
 * @return mat4* 
 */
mat4* rotate_point_x(vec4 point, GLfloat theta);
mat4* rotate_point_x_rad(vec4 point, GLfloat theta);

mat4* rotate_point_y(vec4 point, GLfloat theta);
mat4* rotate_point_y_rad(vec4 point, GLfloat theta);

mat4* rotate_point_z(vec4 point, GLfloat theta);
mat4* rotate_point_z_rad(vec4 point, GLfloat theta);

mat4* rotate_axis(vec4 axis, GLfloat theta);
mat4* rotate_axis_rad(vec4 axis, GLfloat theta);

mat4* rotate_point_axis(vec4 fixed, vec4 axis, GLfloat theta_z);
mat4* rotate_point_axis_rad(vec4 fixed, vec4 axis, GLfloat theta_z);

#endif