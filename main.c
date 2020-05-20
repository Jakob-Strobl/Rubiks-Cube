/*
 * main.c -> Lab06 Executable
 * 
 *  Created on: Oct 31, 2019
 *      Author: Jakob Strobl
 * 
 */

#ifdef __APPLE__    // If on MacOS X device, use these headers..

#define GL_SILENCE_DEPRECATION // Quiet all the warnings because apple stopped supporting OpenGL :(
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>

#else // If on non-Mac OS X operating systems

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>

#endif // __APPLE__

/* Include any required standard libraries here */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include <pthread.h> 

/* Include any custom header files here */
#include "initShader.h" // Thanks Dr. Tan
#include "matrix.h"     // Matrix library -> it points to most updated version in the MAtrix folder
#include "shapes.h"     // Library of various shapes.
// #include "colors.h"     // Library for getting colors 
#include "camera.h"
#include "solve_rc.h"


extern vec4 cube_vertices[36]; // Shapes.c

const int WINDOW_WIDTH = 600;
const int WINDOW_HEIGHT = 600;

// Define bbuffer offset for init()
#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))

// Mouse roataion functions and variables
GLfloat getWorldX(int x) {
    return ((float)x - (WINDOW_WIDTH/2)) / (WINDOW_WIDTH/2);
}

GLfloat getWorldY(int y) {
    return -1 * ((float)y - (WINDOW_HEIGHT/2)) / (WINDOW_HEIGHT/2);
}

GLfloat getWorldZ(GLfloat x, GLfloat y) {
    GLfloat radius = 1;
    GLfloat length = (x*x) + (y*y);
    GLfloat z;

    if (length <= radius) {
        z = sqrtf(radius - length);
    } else {
        z = 0;
    }
    
    return z;
}

struct mouseEvent {
    int button;
    vec4 pos;
    struct mouseEvent *next;
    struct mouseEvent *prev;
} click;


struct mouseEvent *clickStart;
struct mouseEvent *clickEnd;

int isRotating = 0;
vec4 *axis;
GLfloat theta;

const int numCubes = 27;
// Ctm globals
mat4 global_ctm = identity_m4m4_matrix;
mat4 cubie_ctm[numCubes];
GLuint cubie_ctm_location;
GLuint global_ctm_location;

mat4 model_view;
GLuint mv_location;

mat4 projection; 
GLuint projection_location;

// Global user variables
const int numCubies = 1;
const int numVerticesPerCube = 36;
const int numVerticesPerRubiksCube = (6 * (18)) + (3 * 8);
const int num_vertices = (numCubies * numVerticesPerRubiksCube);    // Number of vertices based on number of triangles

// Vertex pointers
vec4 vertices[num_vertices];
vec4 colors[num_vertices];
vec4 normals[num_vertices];

// Rubiks Cube 
//  Colors for a rubiks cube {Green}, {Red},  {Blue},    {Orange},    {White},   {Yellow}
vec4 rubiks_colors[6] = {{0,1,0,1}, {1,0,0,1}, {0,0,1,1}, {1,0.5,0,1}, {1,1,1,1}, {1,1,0,1}};
vec4 edge_color = {0.1,0.1,0.1,1};
//  Positions
//   The values are the index for each cubie_ctm;

vec3 rubiks_pos[numCubes];

// Front goes top down, left to right
// r is the amount of rotation
void rotateFront(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with z positions of 1
        if (rubiks_pos[i][2] == 1) {
            cubie_ctm[i] = *matMult(rotate_z(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_z(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

void rotateRight(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with x positions of 1
        if (rubiks_pos[i][0] == 1) {
            cubie_ctm[i] = *matMult(rotate_x(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_x(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

void rotateBack(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with x positions of 1
        if (rubiks_pos[i][2] == -1) {
            cubie_ctm[i] = *matMult(rotate_z(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_z(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

void rotateLeft(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with x positions of 1
        if (rubiks_pos[i][0] == -1) {
            cubie_ctm[i] = *matMult(rotate_x(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_x(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

void rotateTop(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with x positions of 1
        if (rubiks_pos[i][1] == 1) {
            cubie_ctm[i] = *matMult(rotate_y(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_y(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

void rotateBottom(float r) {
    for (int i = 0; i < numCubes; i++) {
        // Rotate cubes with x positions of 1
        if (rubiks_pos[i][1] == -1) {
            cubie_ctm[i] = *matMult(rotate_y(r), &cubie_ctm[i]);
            // Set new positions of rotated cubes
            matVectorMult_inplace(rotate_y(r), rubiks_pos[i], rubiks_pos[i]);
        }
    }
}

// Animation 
bool isAnimating = 0;
typedef enum rotation {
    NONE = 0,
    FRONT,          // 1
    RIGHT,          // 2
    BACK,           // 3
    LEFT,           // 4
    TOP,            // 5
    BOTTOM,         // 6
    DONE,           // 7    End of shuffling/solution
    ROTATION_NUM    // 8
} face;

face rotating = NONE;
GLfloat rotation_amount = 90;
int rotation_speed = 1;
int rotation_step = 0;
int rotation_total_steps = 50;

bool isShuffling = false;
const int shuffleSize = 30;
face shuffleOrder[shuffleSize];
int shuffleIndex = 0;

bool isGettingSolution = false;
bool isSolving = false;
char *string_solution = NULL;
face *solutionOrder = NULL;
int solutionIndex = 0;

pthread_t thread_id;
void translateSolution();
void *thread_getSolution(void *vargp) {
    string_solution = solve_rc();
            
    //Translate into something I can animate easily
    // printf("Solution: %s\n", string_solution);
    // If the first character is 0, our solution is empty;
    if (string_solution[0] != 0) {
        translateSolution();
        isSolving = true;
    }
    isGettingSolution = false;
}
// The solutionOrder array is terminated by a NONE rotation
void translateSolution() {
    solutionIndex = 0;
    if (string_solution != NULL) {
        int num_turns = 0;
        // Read through string to get size, num at every other character
        for(char *c = string_solution; *c != 0; c += 2) {
            num_turns += *(c + 1) - '0';
        }

        // printf("num turns = %d\n", num_turns);
        int times = 0;
        face rotate = NONE;
        if (num_turns > 0) { 
            solutionOrder = malloc(sizeof(face) * (num_turns + 1));
            for(char *c = string_solution; *c != 0; c += 2) {
                // first (current) character is which face to rotate
                // Second (next) character is how many times to rotate
                times = *(c+1) - '0';
                if(*c == 'F') {
                    rotate = FRONT;
                } else if (*c == 'R') {
                    rotate = RIGHT;
                } else if (*c == 'B') {
                    rotate = BACK;
                } else if (*c == 'L') {
                    rotate = LEFT;
                } else if (*c == 'U') {
                    rotate = TOP;
                } else if (*c == 'D') {
                    rotate = BOTTOM;
                }

                for(int i = 0; i < times; i++, solutionIndex++) {
                    solutionOrder[solutionIndex] = rotate;
                }
            }
            solutionOrder[solutionIndex] = DONE;

            // printf("Translated solution: ");
            // for (int i = 0; ; i++) {
            //     printf("%d ", solutionOrder[i]);
            //     if (solutionOrder[i] == DONE) 
            //         break;
            // }
            // printf("\n");

            solutionIndex = 0;
        }
    }
}

// Init variables
int i = 0;
int j = 0;
int v_index = 0;

// Load the shader files and initialize vertex data for OpenGL to later render.
void init(void) {
    // Initialize the sphere 
    mat4 *temp_ctm = chainCTM(1, scale(0.15, 0.15, 0.15));
    vec4 *cube = cube_vertices;
    cube = getRubiksCube();
    for (i = 0; i < num_vertices; i++) {
        matVectorMult_inplace(temp_ctm, cube[i], vertices[v_index]);
        if (i < 36) {
            colors[v_index][0] = rubiks_colors[i/6][0];
            colors[v_index][1] = rubiks_colors[i/6][1];
            colors[v_index][2] = rubiks_colors[i/6][2];
            colors[v_index][3] = rubiks_colors[i/6][3];
        } else {
            colors[v_index][0] = edge_color[0];
            colors[v_index][1] = edge_color[1];
            colors[v_index][2] = edge_color[2];
            colors[v_index][3] = edge_color[3];
        }
        v_index += 1;
    }

    // Normal per triangle (angle per 3 vertices)
    v_index = 0;
    vec4 v1, v2, normal;
    assert(num_vertices % 3 == 0);
    for (i = 0; i < num_vertices; i+=3) {
        get_vector_inplace(vertices[v_index], vertices[v_index+1], v1); 
        get_vector_inplace(vertices[v_index+1], vertices[v_index+2], v2);
        vecCrossProd_inplace(v1, v2, normal);
        vecNormalize_inplace(normal, normal);
        for (int v = 0; v < 3; v++) {
            normals[v_index][0] = normal[0];
            normals[v_index][1] = normal[1];
            normals[v_index][2] = normal[2];
            normals[v_index][3] = normal[3];
            v_index += 1;
        }
    }


    // Set up the ctms for each cube
    GLfloat x, y, z;
    int x_index, y_index, z_index;
    y_index = 1;
    z_index = 2;
    for (i = 0; i < 27; i++) {
        if (i % 3 == 0) 
            x_index = -1;
        if (i % 3 == 0)
            y_index -= 1;
        if (i % 9 == 0) {
            z_index -= 1;
            y_index = 1;
        }

        x = 0.305 * x_index;
        y = 0.305 * y_index;
        z = 0.305 * z_index;
        // printf("[%d] = %f, %f, %f\n", i, x,  y, z);
        cubie_ctm[i] = *chainCTM(1, translate(x,y,z));
        rubiks_pos[i][0] = x_index;
        rubiks_pos[i][1] = y_index;
        rubiks_pos[i][2] = z_index;
        // printf("POS[%d] = %f, %f, %f\n", i, rubiks_pos[i][0], rubiks_pos[i][1], rubiks_pos[i][2]);
        x_index += 1;
    }

    // Set initial position (rotation) of the rubiks cube
    // global_ctm = *chainCTM(2, rotate_x(30), rotate_y(-35)); 

    r_string_reset();

    free(temp_ctm);

    // Load vertex and frame shader files 
    GLuint program = initShader("vshader.glsl", "fshader.glsl");
    glUseProgram(program);  // Let OpenGL compile shader files

    // For now just use one large array object
    GLuint vao; // Vertex Array Object
    glGenVertexArrays(1, &vao); // 1 vertex array object
    glBindVertexArray(vao);

    // Load arrays you want to render to the screen
    GLuint buffer;
    glGenBuffers(1, &buffer);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors) + sizeof(normals), NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices), sizeof(colors), colors);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(vertices) + sizeof(colors), sizeof(normals), normals);

    // Define attribute vPosition -> defined in .glsl file 
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));

    // Define attribute vColor -> defined in .glsl file 
    GLuint vColor = glGetAttribLocation(program, "vColor");
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(vertices)));

    GLuint vNormal = glGetAttribLocation(program, "vNormal");
    glEnableVertexAttribArray(vNormal);
    glVertexAttribPointer(vNormal, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(vertices) + sizeof(colors)));


    // Initialize location of current transformation matrix
    //  In GLuint program, set a uniform location bound to label 'ctm'
    cubie_ctm_location = glGetUniformLocation(program, "cubie_ctm");
    global_ctm_location = glGetUniformLocation(program, "global_ctm");

    mv_location = glGetUniformLocation(program, "model_view");

    projection_location  = glGetUniformLocation(program, "projection");

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.99, 0.99, 0.99, 1.0); // Set the clear color -> GL_COLOR_BUFFER_BIT
    glDepthRange(1,0);
}

// The function that renders loaded vertices to the screen
// Called each time we want to render a new screen
void display(void) {
    // Clear the current screen in window so we can draw current vertex values.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    glUniformMatrix4fv(global_ctm_location, 1, GL_FALSE, (GLfloat *) &global_ctm);

    glUniformMatrix4fv(mv_location, 1, GL_FALSE, (GLfloat *) &model_view);
    glUniformMatrix4fv(projection_location, 1, GL_FALSE, (GLfloat *) &projection);

    // Draw cubes
    for (i = 0; i < 27; i++) {
        glUniformMatrix4fv(cubie_ctm_location, 1, GL_FALSE, (GLfloat *) &cubie_ctm[i]);
        glDrawArrays(GL_TRIANGLES, 0, num_vertices);  // Draw using triangles, use triangles from index by a number
    }

    // We are using dual buffer rendering
    // Swap old buffer to buffer that just fininished above
    glutSwapBuffers();
}

// Handles user input
void keyboard(unsigned char key, int mousex, int mousey)
{
    if (key == 'q') // If user hits 'q' -> quit program
    	exit(0);

    
    if ( !isAnimating && !isShuffling && !isSolving && !isGettingSolution) {
        if (key == 'f') {
            isAnimating = true;
            rotating = FRONT;
        }
        if (key == 'r') {
            isAnimating = true;
            rotating = RIGHT;
        }
        if (key == 'b') {
            isAnimating = true;
            rotating = BACK;
        }
        if (key == 'l') {
            isAnimating = true;
            rotating = LEFT;
        }
        if (key == 'u') {
            isAnimating = true;
            rotating = TOP;
        }
        if (key == 'd') {
            isAnimating = true;
            rotating = BOTTOM;
        }
        if (key == 's') {
            //shuffle!
            isShuffling = true;
            // Get the shuffle order
            // printf("shuffle order: ");
            for (i = 0; i < shuffleSize - 1; i++) {
                shuffleOrder[i] = (rand() % (ROTATION_NUM - 2)) + 1;
                // printf(" %d", shuffleOrder[i]);
            }
            shuffleOrder[shuffleSize - 1] = DONE; // Terminate with DONE
        }
        // Space key 
        if (key == ' ') {
            //Get the solution!
            isGettingSolution = true;
            pthread_create(&thread_id, NULL, thread_getSolution, NULL);
        }
    }

    // Anytime a user types on keyboard, redraw on the screen.
    glutPostRedisplay(); 
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON: // LMB Click
            // If pressed down 
            if (state == GLUT_DOWN)  {
                isRotating = 0;
                click.button = button;
                click.pos[0] = getWorldX(x); // x
                click.pos[1] = getWorldY(y); // y
                click.pos[2] = getWorldZ(getWorldX(x), getWorldY(y)); // z
                click.pos[3] = 0; // w
                click.next = NULL;
                click.prev = NULL;
                clickStart = &click;
                clickEnd = &click;
            } else { // Mouse button released
                // TODO reset click back to default values.
                // Free the doubly linked list
                if (clickEnd->prev != NULL) {
                    // First Calculate Continual rotation
                    isRotating = 1;
                    for(struct mouseEvent *node = clickEnd->prev; node->prev != NULL; node = node->prev) {
                        free(node->next);
                        node->next = NULL;
                    }
                } else {
                    clickStart = NULL;
                    clickEnd = NULL;
                }
            }
            break;

        case GLUT_RIGHT_BUTTON: // RMB Click
            break;
    }
}

void motion(int x, int y) {
    /*  Track motion  */
    //Allocate memory for next motion 
    clickEnd->next = malloc(sizeof(struct mouseEvent));

    // Set next new node prev to the current end of the list. 
    clickEnd->next->prev = clickEnd;

    // Linked List end is now point to next node. 
    clickEnd = clickEnd->next;

    // Inherit button value from the first node (start) in the list 
    clickEnd->button = clickStart->button;

    // Set position
    clickEnd->pos[0] = getWorldX(x);
    clickEnd->pos[1] = getWorldY(y);
    clickEnd->pos[2] = getWorldZ(getWorldX(x), getWorldY(y));
    clickEnd->pos[3] = 0;

    // printf("\n\nCurrent Position: ");
    // vecPrint(clickEnd->pos);

    // Set next to null
    clickEnd->next = NULL;

    // Get the axis we want to rotate about
    axis = vecNormalize(*vecCrossProd(clickEnd->prev->pos, clickEnd->pos));
    
    // Get dot product to calculate theta
    GLfloat dp = vecDotProd(clickEnd->prev->pos, clickEnd->pos);

    if (dp > 1)
        theta = acosf(1);
    else if (dp < -1) 
        theta = acosf(-1);
    else
        theta = acosf(dp);


    global_ctm = *matMult(rotate_axis_rad(*axis, theta), &global_ctm);

}

void idle(void) {
    // Rotate about a arbitrary axis
    if (isRotating && fabs(theta * 180 / M_PI) > 1) {
        global_ctm = *matMult(rotate_axis_rad(*axis, theta), &global_ctm);
    } else {
        isRotating = 0;
    }

    // If we are shuffling lets go through the randomly generated turns
    if (isShuffling) {
        if ( !isAnimating ) {
            #ifdef __APPLE__
            rotation_total_steps = 10;
            #else
            rotation_total_steps = 20;
            #endif
            rotating = shuffleOrder[shuffleIndex];
            shuffleIndex += 1;
            isAnimating = true;
        }
    } else if (isSolving) {
        if ( !isAnimating ) {
            #ifdef __APPLE__
            rotation_total_steps = 10;
            #else
            rotation_total_steps = 20;
            #endif
            rotating = solutionOrder[solutionIndex];
            solutionIndex += 1;
            isAnimating = true;
        }
    }

    // Process rotation animations 
    if (isAnimating) {
        float r = rotation_amount * rotation_speed / rotation_total_steps;
        if (rotating == FRONT) {
            if (rotation_step == 0)
                r_string_front();

            rotateFront(-r);
        } else if (rotating == RIGHT) {
            if (rotation_step == 0)
                r_string_right();

            rotateRight(-r);
        } else if (rotating == BACK) {
            if (rotation_step == 0)
                r_string_back();
             
            rotateBack(r);
        } else if (rotating == LEFT) {
            if (rotation_step == 0)
                r_string_left();
             
            rotateLeft(r);
        } else if (rotating == TOP) {
            if (rotation_step == 0)
                r_string_up();
             
            rotateTop(-r);
        } else if (rotating == BOTTOM) {
            if (rotation_step == 0)
                r_string_down();
             
            rotateBottom(r);
        } else if (rotating == NONE) {
            printf("WARNING: Is animating without a direction :(\n");
        }

        rotation_step += rotation_speed;

        if (rotating == DONE) {
            isAnimating = false;
            isSolving = false;
            isShuffling = false;
            #ifdef __APPLE__
            rotation_total_steps = 25;
            #else
            rotation_total_steps = 40;
            #endif
            shuffleIndex = 0;
            rotation_step = 0;
            solutionIndex = 0;
            if (solutionOrder != NULL) {
                free(solutionOrder);
                solutionOrder = NULL;
            }
        }
    }

    // Once rotation animation is complete, set isAnimating to false;
    if (rotation_step == rotation_total_steps) {
        isAnimating = false;
        rotating = NONE;
        rotation_step = 0;
        for(i = 0; i < numCubes; i++) {
            // Not a perfect idea, but its simple and it works
            // round the position at end of an rotation, so we don't have to worry about floating point error
            rubiks_pos[i][0] = roundf(rubiks_pos[i][0]);
            rubiks_pos[i][1] = roundf(rubiks_pos[i][1]);
            rubiks_pos[i][2] = roundf(rubiks_pos[i][2]);
        }
    } else if (rotation_amount > 90) {
        printf("Warning: Rotation overshot!!!!\n");
    }

    // Ask OpenGL to redraw the window.
    glutPostRedisplay(); 
}

// Main process to initialize and set up 
int main(int argc, char **argv) {
    // Seed the random generator to current system time.
    srand(time(NULL));

    glutInit(&argc, argv);

    // Set color mode RGBA, Set as double render buffer, Enable depth buffer
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH); 
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);       // Set window size
    glutInitWindowPosition(100,100);    // Set position
    glutCreateWindow("Rubik's Cube");  // Set window name
    
    model_view = *look_at((vec4){0, 0, 1.5, 0}, (vec4){0, 0, 0, 0}, (vec4){0, 1, 0, 0});
    projection = *perspective(90, 1, -0.1, -9);
    
    #ifndef __APPLE__   // If using macOS X, skip glewInit()
        glewInit();
    #endif  // __APPLE__
    
    init(); // Call user defined function above 

    // Set the function to call when we want to redraw the screen.
    glutDisplayFunc(display);  // display signature = void display(void);

    // Set function to invoke on keyboard event
    glutKeyboardFunc(keyboard);

    // Set mouse event function
    glutMouseFunc(mouse);

    // Set motion event function
    glutMotionFunc(motion);

    // Set function to invoke when idle
    glutIdleFunc(idle);

    // Start the main render loop of OpenGL
    glutMainLoop();

    return 0;
}