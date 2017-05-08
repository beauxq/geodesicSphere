#include <iostream>

#include "gs.h"

int main(int argc, char* argv[]) {
    gs::VertexList vertices;
    gs::TriangleList triangles;

    gs::geodesic_sphere(vertices, triangles, 1, 5);

    gs::create_stl("my15.stl", vertices, triangles);

    return 0;
}