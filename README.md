# Rubiks-Cube

 Rubik's Cube Written in C with OpenGL from scratch.

## About

This was the final project for my Computer Graphics class during my Computer Science Undergrad at Pitt. Majority of the source code was written by me, the rest was provided by my professor (see below). This project combined all the topics we learned during the course: Generating and loading geometry, animation (tied to FPS :/ - time constraints), arcball rotation with the mouse, perspective projection, surface normals, Phong's lighting model / some minor shader programming.

## [Demo: Video Link](https://www.youtube.com/watch?v=98dx05yLDUw)

## Controls

The controls are based on the intial perspective of the rubiks cube when you start the program. 
So the front face is based on the angle when you first start the program, not the face you see after rotating.

### Q - Quit Program

### F - Move Front Face (Green side)

### D - Rotate Bottom Face

### B - Rotate Back Face

### U - Rotate Top Face

### R - Rotate Right Face

### L - Rotate Left Face

### S - Shuffles the cube

### ' ' (space) - Solves the cube

### Clink and Drag with the Mouse - Rotate cube

- Click, drag, release for endless spinngin :)

## Source Code (Written By Me)

readme: this file
makefile: make file
camera.c/.h: Code for orthographic and perspective functions.
matrix.c/.h: All matrix operations and definitions
shapes.c/.h: Functions that make 3D shapes
main.c: The file that contains main().
fshader.glsl: Frame shader for textures.
vshader.glsl: Vertex shader for ctm, model view, and perspective matrix.

### Source Code Not Written by Me (Files given by Dr. Tan)

- initShader.c
  - Functions to read, compile, and load shader (.glsl) files.
- initShader.h
  - Header file for initShader.c
- solve_rc.c/solve_rc.h
  - Rubiks cube solver adapted to C by Dr. Tan

## Running The Program Locally

To run the lab:

- If you are using mac, type
  - make mac (You'll need the necessary dependencies/config of library files D: )
  - ./main (start program)
- If you are using linux, type
  - make linux
  - ./main (start program)
  - I did not test with linux for this project. I didn't have enough time for extensive testing.