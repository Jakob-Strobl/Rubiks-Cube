/*
 * shapes.c 
 *      This file contains methods to create a variety of shapes.
 *          - Circle    (lab03)
 *          - Cone      (lab03)
 * 
 *  Created on: Sep 17, 2019
 *      Author: Jakob Strobl
 * 
 */

#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "shapes.h"

vec4 cube_vertices[36] = {{-1, -1, 1, 1}, {1, -1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {-1, 1, 1, 1}, {-1, -1, 1, 1}, // Front
                 {1, -1, -1, 1}, {-1, -1, -1, 1}, {-1, 1, -1, 1}, {-1, 1, -1, 1}, {1, 1, -1, 1}, {1, -1, -1, 1},    // Back
                 {1, -1, 1, 1}, {1, -1, -1, 1}, {1, 1, -1, 1}, {1, 1, -1, 1}, {1, 1, 1, 1}, {1, -1, 1, 1},          // Right
                 {-1, -1, -1, 1}, {-1, -1, 1, 1}, {-1, 1, 1, 1}, {-1, 1, 1, 1}, {-1, 1, -1, 1}, {-1, -1, -1, 1},    // Left
                 {-1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, -1, 1}, {1, 1, -1, 1}, {-1, 1, -1, 1}, {-1, 1, 1, 1},          // Top
                 {1, -1, 1, 1}, {-1, -1, 1, 1}, {-1, -1, -1, 1}, {-1, -1, -1, 1}, {1, -1, -1, 1}, {1, -1, 1, 1}};   // Bot

typedef struct {
    triangle triBot;
    triangle triTop;
} face;


/**
 * @brief Get a Rectangle based of the parameters
 * 
 * @param pointA - x, y, z position of the bottom-left corner.
 * @param pointB - x, y, z position of the bottom-right corner.
 * @param pointC - x, y, z position of the top-right corner.
 * @param pointD - x, y, z position of the top-left corner.
 * @return vec4* - Vertex array of 2 triangles (6 vertices) to form the rectangle
 */
vec4* getRectangle(vec4 pointA, vec4 pointB, vec4 pointC, vec4 pointD) {
    vec4 *rectangle = malloc(sizeof(triangle) * 2); // Two triangles per rectangle
    triangle *tri = (triangle *)rectangle;

    // First Triangle
    // Set point A of triangle
    tri->A[0] = pointA[0];
    tri->A[1] = pointA[1];
    tri->A[2] = pointA[2];
    tri->A[3] = 1;

    // Set point B
    tri->B[0] = pointB[0];
    tri->B[1] = pointB[1];
    tri->B[2] = pointB[2];
    tri->B[3] = 1;

    // Set point C
    tri->C[0] = pointC[0];
    tri->C[1] = pointC[1];
    tri->C[2] = pointC[2];
    tri->C[3] = 1;

    // printf("Triangle A: \n");
    // vecPrint(tri->A);
    // vecPrint(tri->B);
    // vecPrint(tri->C);

    // Next triangle!
    tri += 1;

    // Set point A
    tri->A[0] = pointC[0];
    tri->A[1] = pointC[1];
    tri->A[2] = pointC[2];
    tri->A[3] = 1;

    // Set point B
    tri->B[0] = pointD[0];
    tri->B[1] = pointD[1];
    tri->B[2] = pointD[2];
    tri->B[3] = 1;

    // Set point C
    tri->C[0] = pointA[0];
    tri->C[1] = pointA[1];
    tri->C[2] = pointA[2];
    tri->C[3] = 1;

    // printf("Triangle B: \n");
    // vecPrint(tri->A);
    // vecPrint(tri->B);
    // vecPrint(tri->C);

    return rectangle;
}

vec4* createPrismWackyPointers(GLfloat length, GLfloat width, GLfloat depth) {
    GLfloat l = length/2;
    GLfloat w = width/2;
    GLfloat d = depth/2;
    vec4 *prism = (vec4*)malloc(sizeof(vec4) * 36);
    face *prism_face = (face *)prism;
    
    // Front Face
    vec4 *front;
    front = getRectangle((vec4){-w, -l, d, 1},(vec4){w, -l, d, 1},(vec4){w, l, d, 1},(vec4){-w, l, d, 1});
    *prism_face = *(face *)front;
    free(front);
    prism_face += 1;

    // Right Face
    vec4 *right;
    right = getRectangle((vec4){w, -l, d, 1},(vec4){w, -l, -d, 1},(vec4){w, l, -d, 1},(vec4){w, l, d, 1});
    *prism_face = *(face *)right;
    free(right);
    prism_face += 1;

    // Back Face
    vec4 *back;
    back = getRectangle((vec4){w, -l, -d, 1},(vec4){-w, -l, -d, 1},(vec4){-w, l, -d, 1},(vec4){w, l, -d, 1});
    *prism_face = *(face *)back;
    free(back);
    prism_face += 1;

    // Left Face
    vec4 *left;
    left = getRectangle((vec4){-w, -l, -d, 1},(vec4){-w, -l, d, 1},(vec4){-w, l, d, 1},(vec4){-w, l, -d, 1});
    *prism_face = *(face *)left;
    free(left);
    prism_face += 1;

    // Top Face
    vec4 *top;
    top = getRectangle((vec4){-w, l, d, 1},(vec4){w, l, d, 1},(vec4){w, l, -d, 1},(vec4){-w, l, -d, 1});
    *prism_face = *(face *)top;
    free(top);
    prism_face += 1;

    // Bottom Face
    vec4 *bot;
    bot = getRectangle((vec4){w, -l, d, 1},(vec4){-w, -l, d, 1},(vec4){-w, -l, -d, 1},(vec4){w, -l, -d, 1});
    *prism_face = *(face *)bot;
    free(bot);
    prism_face += 1;

    return prism;
}

vec4* createPrismMath(GLfloat length, GLfloat width, GLfloat depth) {
    // Two triangles per face of a prism
    // 6 faces in a rectangular prism
    vec4 *prism = malloc(sizeof(face) * 6);
    vec4 *prism_ptr = prism;
    vec4 pointA, pointB, pointC, pointD;

    pointA[0] = -width/2;
    pointA[1] = -length/2;
    pointA[2] = depth/2;
    pointA[3] = 1;

    pointB[0] = width/2;
    pointB[1] = -length/2;
    pointB[2] = depth/2;
    pointB[3] = 1;

    pointC[0] = width/2;
    pointC[1] = length/2;
    pointC[2] = depth/2;
    pointC[3] = 1;

    pointD[0] = -width/2;
    pointD[1] = length/2;
    pointD[2] = depth/2;
    pointD[3] = 1;

    vec4 *front = getRectangle(pointA, pointB, pointC, pointD);
    vec4 *front_ptr = (vec4 *)front;

    // Rotate faces to front, right, left, bottom
    for (int f = 0; f < 4; f++) {
        for (int i = 0; i < 6; i++) {
            matVectorMult_inplace(rotate_y(90*(f)), *front_ptr, (GLfloat *)prism_ptr);
            vecPrintAll(*prism_ptr);
            front_ptr += 1;
            prism_ptr += 1;
        }
        front_ptr = front;
    }
    
    // Rotate face to bottom
    front_ptr = front;
    for (int i = 0; i < 6; i++) {
        matVectorMult_inplace(rotate_x(90), *front_ptr, (GLfloat *)prism_ptr);
        front_ptr += 1;
        prism_ptr += 1;
    }

    // Rotate face to top
    front_ptr = front;
    for (int i = 0; i < 6; i++) {
        matVectorMult_inplace(rotate_x(-90), *front_ptr, (GLfloat *)prism_ptr);
        front_ptr += 1;
        prism_ptr += 1;
    }

    free(front);
    return prism;
}


vec4* getRubiksCube() {
    // 6 faces and 12 sides -> each are a rectangle
    vec4 *cube = malloc(((sizeof(vec4) * 6) * (18))  + ((sizeof(vec4) * 3) * 8));
    int cubeIndex = 0;
    vec4 *face, *edge;
    GLfloat face_w = 0.9;
    GLfloat face_h = 0.9;
    GLfloat gap = 0.1;
    GLfloat edge_w = 0.1;
    GLfloat edge_h = 0.9;

    // Get the front face
    face = getRectangle((vec4){-face_w, -face_h, face_w + gap, 1},
                        (vec4){ face_w, -face_h, face_w + gap, 1},
                        (vec4){ face_w,  face_h, face_w + gap, 1},
                        (vec4){-face_w,  face_h, face_w + gap, 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Get the right face
    face = getRectangle((vec4){ face_w + gap, -face_h,  face_w, 1},
                        (vec4){ face_w + gap, -face_h, -face_w, 1},
                        (vec4){ face_w + gap,  face_h, -face_w, 1},
                        (vec4){ face_w + gap,  face_h,  face_w, 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Get the back face
    face = getRectangle((vec4){ face_w, -face_h, -(face_w + gap), 1},
                        (vec4){-face_w, -face_h, -(face_w + gap), 1},
                        (vec4){-face_w,  face_h, -(face_w + gap), 1},
                        (vec4){ face_w,  face_h, -(face_w + gap), 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Get the left face
    face = getRectangle((vec4){ -(face_w + gap), -face_h, -face_w, 1},
                        (vec4){ -(face_w + gap), -face_h,  face_w, 1},
                        (vec4){ -(face_w + gap),  face_h,  face_w, 1},
                        (vec4){ -(face_w + gap),  face_h, -face_w, 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Get the top face
    face = getRectangle((vec4){ -face_w, (face_h + gap),  face_w, 1},
                        (vec4){  face_w, (face_h + gap),  face_w, 1},
                        (vec4){  face_w, (face_h + gap), -face_w, 1},
                        (vec4){ -face_w, (face_h + gap), -face_w, 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Get the bottom face
    face = getRectangle((vec4){  face_w, -(face_h + gap),  face_w, 1},
                        (vec4){ -face_w, -(face_h + gap),  face_w, 1},
                        (vec4){ -face_w, -(face_h + gap), -face_w, 1},
                        (vec4){  face_w, -(face_h + gap), -face_w, 1});
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = face[i][j];
        }
        cubeIndex += 1;
    }
    free(face);

    // Start the edges
    // cube[1] is bott right of front face, cube[6] is bot left of right face, etc...
    // Edge betweeen front and right face
    edge = getRectangle(cube[1], cube[6], cube[10], cube[2]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between right and back face
    edge = getRectangle(cube[7], cube[12], cube[16], cube[8]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between back and left face
    edge = getRectangle(cube[13], cube[18], cube[22], cube[14]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between left and right 
    edge = getRectangle(cube[19], cube[0], cube[4], cube[20]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between front and top
    edge = getRectangle(cube[4], cube[3], cube[25], cube[24]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between right and top
    edge = getRectangle(cube[10], cube[9], cube[26], cube[25]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between back and top
    edge = getRectangle(cube[16], cube[15], cube[28], cube[26]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between left and top
    edge = getRectangle(cube[22], cube[21], cube[29], cube[28]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between front and bot
    edge = getRectangle(cube[1], cube[0], cube[31], cube[30]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between right and bot
    edge = getRectangle(cube[7], cube[6], cube[30], cube[34]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Edge between back and bot
    edge = getRectangle(cube[13], cube[12], cube[34], cube[32]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);


    // Edge between left and bot
    edge = getRectangle(cube[19], cube[18], cube[32], cube[31]);
    for(int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            cube[cubeIndex][j] = edge[i][j];
        }
        cubeIndex += 1;
    }
    free(edge);

    // Fill in the triangle gaps
    // Top Front-Right
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[3][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[10][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[25][j];
    }
    cubeIndex += 1;

    // Top Back-Right
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[9][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[16][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[26][j];
    }
    cubeIndex += 1;

    // Top Back-Left
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[15][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[22][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[28][j];
    }
    cubeIndex += 1;

    // Top Front-Left
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[21][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[4][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[29][j];
    }
    cubeIndex += 1;

    // Bot Front-Right
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[1][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[30][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[6][j];
    }
    cubeIndex += 1;

    // Bot Back-Right
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[7][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[34][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[12][j];
    }
    cubeIndex += 1;

    // Bot Back-left
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[13][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[33][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[18][j];
    }
    cubeIndex += 1;

    // Bot Front-left
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[19][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[31][j];
    }
    cubeIndex += 1;
    for (int j = 0; j < 4; j++) {
        cube[cubeIndex][j] = cube[0][j];
    }
    cubeIndex += 1;

    return cube;
}

vec4* getSphere(int numRings, int numSides) {
    size_t rectangle = sizeof(triangle) * 2;
    vec4 *sphere = malloc(rectangle * numRings * numSides);
    int sphereIndex = 0;

    // Thinnk of a 2D grid of longitude and latitude
    //  +1 is so we wrap around 360 degrees
    vec4 map[numRings+1][numSides+1];

    // Radius of circle is always 1 (fill i)
    GLfloat radius = 1;

    // Theta is used for the angle for logitude (rings)
    // The azimuth
    GLfloat theta = 2 * M_PI / numRings;

    // Phi is used for the angle for latitude (sides)
    // The inclination
    GLfloat phi = M_PI / numSides;

    // Calculate the x, y, z points of the sphere 
    for (int ring = 0; ring < numRings+1; ring++) {
        GLfloat longitude = ring * theta;
        for (int side = 0; side < numSides+1; side++) {
            GLfloat latitude = side * phi;
            GLfloat x = radius * sinf(latitude) * cosf(longitude);
            GLfloat y = radius * sinf(latitude) * sinf(longitude);
            GLfloat z = radius * cosf(latitude);
            map[ring][side][0] = x;
            map[ring][side][1] = y;
            map[ring][side][2] = z;
            map[ring][side][3] = 1;
        }
    }


    // Figure out the triangles from the 2D map of points
    for (int ring = 0; ring < numRings; ring++) {
        for (int side = 0; side < numSides; side++) {
            // Draw the rectangle! (using two triangles)
            
            // Triangle A
            // Top left
            for (int p = 0; p < 4; p++) { 
                sphere[sphereIndex][p] = map[ring][side][p];
            }
            sphereIndex+=1;
            // Bot left
            for (int p = 0; p < 4; p++) {
                sphere[sphereIndex][p] = map[ring][side+1][p];
            }
            sphereIndex+=1;
            // Bot right
            for (int p = 0; p < 4; p++) {
                sphere[sphereIndex][p] = map[ring+1][side+1][p];
            }
            sphereIndex+=1;

            // Triangle B
            // Bot right
            for (int p = 0; p < 4; p++) {
                sphere[sphereIndex][p] = map[ring+1][side+1][p];
            }
            sphereIndex+=1;
            // Top right
            for (int p = 0; p < 4; p++) {
                sphere[sphereIndex][p] = map[ring+1][side][p];
            }
            sphereIndex+=1;
            // Top left
            for (int p = 0; p < 4; p++) {
                sphere[sphereIndex][p] = map[ring][side][p];
            }
            sphereIndex+=1;
        }
    }

    return sphere;
}

// Project 1
// Returns # of vertices = numRings * numSides (numRectangles) * 3 * 2
vec4* renderTorus(GLfloat R, GLfloat r, int numRings, int numSides) {
    size_t rectangle = sizeof(triangle) * 2;
    vec4 *torus = malloc(rectangle * numSides * numRings);
    triangle *torus_ptr = (triangle *)torus;

    GLfloat x, y, z, s, t, h;
    GLfloat theta = M_PI * 2;

    for (int ring = 0; ring < numRings; ring++) {
        for (int side = 0; side < numSides; side++) {
            s = (ring + 0) % numRings + (R / 2);
            t = (side + 0) % numSides;
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / numSides);
            y = r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / numSides);
            vec4 pointA = {x, y, z, 1}; //wtf - + z

            s = (ring + 1) % numRings + (R / 2);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / numSides);
            y = r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / numSides);
            vec4 pointB = {x, y, z, 1};

            s = (ring + 1) % numRings + (R / 2);
            t = (side + 1) % numSides;
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / numSides);
            y = r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / numSides);
            vec4 pointC = {x, y, z, 1}; //wtf - + z

            s = (ring + 0) % numRings + (R / 2);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / numSides);
            y = r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / numSides);
            vec4 pointD = {x, y, z, 1};

            vec4 *rect_vertices = getRectangle(pointA, pointB, pointC, pointD); // return 6 vertices
            triangle *rect_tri = (triangle *)rect_vertices;
            *torus_ptr = *rect_tri; // Copy first triangle from rect.

            rect_tri += 1;
            torus_ptr += 1;

            *torus_ptr = *rect_tri; // Copy second triangle from rect.

            // Go to next triangle for vertices
            torus_ptr += 1;
            free(rect_vertices);
        }
    }
    return torus;
}

vec4* renderSpring(GLfloat R, GLfloat r, int numRings, int numSides, int numRotations) {
    printf("Sides %d   |  Rings %d", numSides, numRings);
    size_t rectangle = sizeof(triangle) * 2;
    vec4 *torus = malloc(rectangle * (numSides+ 2)* numRings);
    triangle *torus_ptr = (triangle *)torus;

    GLfloat x, y, z, s, t, h;
    GLfloat theta = M_PI * 2;

    h = 0.1 * numRotations;

    for (int ring = 0; ring < numRings; ring++) {
        
        if (ring == 0) {
            // Make circle ends to the spring
            GLfloat *ptr = (GLfloat *)torus;

            // Calculate each triangle's 3 points (vec4)
            for (int i = 0; i < numRings; i++) {
                s = (i + 0) + (R / 2);
                t = 0;

                /* Calculate point A (first) of triangle (TOP) */
                // Calculate x position
                *ptr = R + r * cosf(s * theta / numRings);
                ptr += 1;  
                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
                ptr += 1;   
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = 0;
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;
                ptr += 1;  

                /* Calculate point B (second) of traingle (LEFT) */
                // Calculate x position
                *ptr = R;
                ptr += 1;  
                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides);
                ptr += 1;   
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = 0;
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;
                ptr += 1;   

                s = (i + 1) + (R / 2);
                /* Calculate point C (third) of triangle (RIGHT) */
                // Calculate x position
                *ptr = R + r * cosf(s * theta / numRings);
                ptr += 1;  
                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
                ptr += 1;   
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = 0;
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;
                ptr += 1;  
            }
            torus_ptr = (triangle *)ptr;
        }
        for (int side = 0; side < numSides; side++) {
            
            // Calculate point A of rectangle 
            s = (ring + 0) + (R / 2);
            t = (side + 0);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
            y = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
            vec4 pointA = {x, y, z, 1}; //wtf - + z

            // Calculate point B of rectangle 
            s = (ring + 1) + (R / 2);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
            y = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
            vec4 pointB = {x, y, z, 1};

            // Calculate point C of rectangle 
            s = (ring + 1) + (R / 2);
            t = (side + 1);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
            y = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
            vec4 pointC = {x, y, z, 1}; //wtf - + z

            // Calculate point D of rectangle 
            s = (ring + 0) + (R / 2);
            x = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
            y = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
            z = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
            vec4 pointD = {x, y, z, 1};

            // Get rectangle and copy over the two triangles to our spring array
            vec4 *rect_vertices = getRectangle(pointA, pointB, pointC, pointD); // return 6 vertices
            triangle *rect_tri = (triangle *)rect_vertices;
            *torus_ptr = *rect_tri; // Copy first triangle from rect.

            rect_tri += 1;
            torus_ptr += 1;

            *torus_ptr = *rect_tri; // Copy second triangle from rect.
            // Go to next triangle for vertices
            torus_ptr += 1;
            free(rect_vertices);
        }
        if (ring == numRings - 1) {
            // Make circle ends to the spring
            GLfloat *ptr = (GLfloat *)torus_ptr;
            int side = numSides-1;
            GLfloat center_z = (R + r * cosf(s * theta / numRings / 4)) * sinf(t * theta / (numSides / numRotations));
            printf("Center %f",center_z );
            // Calculate each triangle's 3 points (vec4)
            for (int i = 0; i < numRings; i++) {
                t = (side + 1);
                s = (i + 0) + (R / 2);

                /* Calculate point A (first) of triangle (TOP) */ 
                // Calculate x position
                *ptr = R * cosf(t * theta / (numSides / numRotations));
                ptr += 1;  
                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides);
                ptr += 1;   
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = center_z;
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;
                ptr -= 3;
                // printf("Point A: ");
                // vecPrint((vec4*)ptr);
                ptr += 3;
                ptr += 1;  

                s = (i + 0) + (R / 2);

                /* Calculate point B (second) of traingle (LEFT) */
                // Calculate x position
                *ptr = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
                ptr += 1;  
                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
                ptr += 1; 
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;

                ptr -= 3;
                // printf("Point B: ");
                // vecPrint((vec4*)ptr);
                ptr += 3;
                ptr += 1;  

                s = (i + 1) + (R / 2);
                /* Calculate point C (third) of triangle (RIGHT) */
                // Calculate x position
                *ptr = (R + r * cosf(s * theta / numRings)) * cosf(t * theta / (numSides / numRotations));
                ptr += 1;  

                // Calculate y position   
                *ptr = h * ((t - numSides/2) * theta / numSides) + r * sinf(s * theta / numRings);
                ptr += 1;   
                // Calculate z position, cirlces are 2D, so Z poisition is zero
                *ptr = (R + r * cosf(s * theta / numRings)) * sinf(t * theta / (numSides / numRotations));
                ptr += 1; 
                // Set w = 1, since this is a line
                *ptr = 1;

                ptr -= 3;
                // printf("Point C: ");
                // vecPrint((vec4*)ptr);
                ptr += 3;
                
                ptr += 1;  
            }
        }
    }
    return torus;
}

/** @funtion renderCircle
 *      Renders a cirlce around the origin using number of triangles defined by parameter numTriangles 
 * 
 *  @returns vec4 *
 *      Returns an array of vertices of triangles
 *      array size = vec4[3 * numTriangles]
 */
vec4* renderCircle(int numTriangles) {
    GLfloat radius = 1;
    // Each triangle needs 3 4x1 vectors, so we multiply number of triangles by 3.
    vec4 *vertices = (vec4*)malloc(sizeof(vec4) * numTriangles * 3);

    // Calculate theta for the angle of each triangle
    GLfloat theta = 2 * M_PI / numTriangles;

    // Start pointer at beginning of array
    GLfloat *ptr = vertices[0];

    // Calculate each triangle's 3 points (vec4)
    for (int i = 0; i < numTriangles; i++) {
        /* Calculate point A (first) of triangle (TOP) */
        // Calculate x position
        *ptr = 0;
        ptr += 1;  
        // Calculate y position   
        *ptr = 0;
        ptr += 1;   
        // Calculate z position, cirlces are 2D, so Z poisition is zero
        *ptr = 0;
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;   

        /* Calculate point B (second) of traingle (LEFT) */
        // Calculate x position
        *ptr = radius * cos(i * theta);
        ptr += 1;  
        // Calculate y position   
        *ptr = radius * sin(i * theta);
        ptr += 1;   
        // Calculate z position, cirlces are 2D, so Z poisition is zero
        *ptr = 0;
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;  

        /* Calculate point C (third) of triangle (RIGHT) */
        // Calculate x position
        *ptr = radius * cos((i + 1) * theta);
        ptr += 1;  
        // Calculate y position   
        *ptr = radius * sin((i + 1) * theta);
        ptr += 1;   
        // Calculate z position, cirlces are 2D, so Z poisition is zero
        *ptr = 0;
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;  
    }

    return vertices;
}

/** @funtion renderCone
 *      Calculate a 3D cone, where the base sits on the x-z plane and points up parallel to y axis.
 *      The origin is at the center of the cone.
 *      You can also modify the size of the cone, by changing parameter size.
 * 
 *  @returns vec4 *
 *      Returns an array of vertices of triangles
 *      array size = vec4[3 * numTriangles]
 */
vec4* renderCone(int numTriangles) {
    // Each triangle needs 3 4x1 vectors, so we multiply number of triangles by 3.
    // Since a cone has two sets
    vec4 *vertices = malloc(sizeof(vec4) * numTriangles * 3 * 2);

    // Calculate theta for the angle of each triangle
    GLfloat theta = 2 * M_PI / (numTriangles / 2);
    //printf("Theta: %f", theta);

    // Start pointer at beginning of array
    GLfloat *ptr = vertices[0];

    // Create the base circle
    for (int i = 0; i < numTriangles/2; i++) {
        /* Calculate point A (first) of triangle (TOP) */
        // X position
        *ptr = 0;
        ptr += 1;  
        // Y position   
        *ptr = -1;
        ptr += 1;   
        // Z position
        *ptr = 0;
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;   

        /* Calculate point B (second) of traingle */
        // X position
        *ptr = cos(i * theta);
        ptr += 1;  
        // Y position   
        *ptr = -1;
        ptr += 1;   
        // Z position
        *ptr = sin(i * theta);;
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;  

        /* Calculate point C (third) of triangle */
        // X position
        *ptr = cos((i + 1) * theta);
        ptr += 1;  
        // Y position   
        *ptr = -1;
        ptr += 1;   
        // Z position
        *ptr = sin((i + 1) * theta);
        ptr += 1; 
        // Set w = 1, since this is a line
        *ptr = 1;
        ptr += 1;  
    }

    // Create the Cone (?) part
    for (int i = 0; i < numTriangles/2; i++) {
        /* Calculate point A of triangle 
         *  TOP: 
         *      x = 0, 
         *      y = height/2, 
         *      z = 0
         */
        
        *ptr = 0;           // X position
        ptr += 1;  

        *ptr = 1;           // Y position   
        ptr += 1;   

        *ptr = 0;           // Z position
        ptr += 1; 
        
        *ptr = 1;
        ptr += 1;  

        /* Calculate point B of triangle 
         *  BOTTOM LEFT: 
         *      x = height/2 * cos(phi + theta), 
         *      y = -1 * height/2, 
         *      z = height/2 * sin(phi + theta)
         */
        
        *ptr = cos((i + 1) * theta); // X position
        ptr += 1;  

        *ptr = -1;                   // Y position 
        ptr += 1;   
        
        *ptr = sin((i + 1) * theta); // Z position
        ptr += 1; 

        *ptr = 1;   // Set w = 1
        ptr += 1;  

        /* Calculate point C (first) of triangle 
         *  BOTTOM RIGHT: 
         *      x = height/2 * cos(phi), 
         *      y = -1 * height/2, 
         *      z = height/2 * sin(phi)
         */
        
        *ptr = cos(i * theta);;
        ptr += 1;  
           
        *ptr = -1;
        ptr += 1;   
        
        *ptr = sin(i * theta);;
        ptr += 1; 
        
        *ptr = 1;
        ptr += 1;   
    }

    return vertices;

}