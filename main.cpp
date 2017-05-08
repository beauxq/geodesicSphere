#include <iostream>

#include "gs.h"

int main(int argc, char* argv[]) {
    gs::VertexList vertices;
    gs::TriangleList triangles;

    gs::geodesic_sphere(vertices, triangles, 1, 5);

    gs::create_binary_stl("mybinary15.stl", vertices, triangles);

    gs::geodesic_sphere(vertices, triangles, 1, 1);

    gs::create_ascii_stl("myascii11stl", vertices, triangles);

    return 0;
}