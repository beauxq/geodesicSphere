// library to create geodesic sphere

#ifndef GEODESICSPHERE_GS_H
#define GEODESICSPHERE_GS_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>

namespace gs {
    const size_t ORIGINAL_VERTEX_COUNT = 12;
    const size_t ORIGINAL_TRIANGLE_COUNT = 20;

    typedef std::vector< double > Vertex;
    typedef std::vector < size_t > Triangle;

    typedef std::vector< Vertex > VertexList;
    typedef std::vector< Triangle > TriangleList;

    Vertex operator+ (const Vertex& lhs, const Vertex& rhs) {
        return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
    }

    Vertex& operator+= (Vertex& lhs, const Vertex& rhs) {
        lhs = {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
        return lhs;
    }

    Vertex operator- (const Vertex& lhs, const Vertex& rhs) {
        return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
    }

    Vertex operator* (const Vertex& lhs, const double& rhs) {
        return {lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs};
    }

    Vertex operator/ (const Vertex& lhs, const int& rhs) {
        return {lhs[0] / (double)rhs, lhs[1] / (double)rhs, lhs[2] / (double)rhs};
    }

    double distance_between_2_vertices(const Vertex& a, const Vertex& b) {
        Vertex diff = a - b;
        return sqrt(diff[0] * diff[0] +
                    diff[1] * diff[1] +
                    diff[2] * diff[2]);
    }

    /**
     * all the vertices of an icosahedron
     * with y axis being up, looking down,
     * starts at top, goes clockwise around next level,
     * then clockwise around next level, then bottom
     */
    const VertexList ORIGINAL_VERTICES = {
            {0.000000, 1.000000, 0.000000},
            {0.894450, 0.447168, 0.000000},
            {0.276400, 0.447168, 0.850673},
            {-0.723625, 0.447168, 0.525744},
            {-0.723625, 0.447168, -0.525745},
            {0.276400, 0.447168, -0.850673},
            {0.723625, -0.447168, 0.525745},
            {-0.276400, -0.447168, 0.850673},
            {-0.894450, -0.447168, -0.000000},
            {-0.276400, -0.447168, -0.850673},
            {0.723625, -0.447168, -0.525745},
            {0.000000, -1.000000, 0.000000}
    };

    bool testDistanceOfAllVertices() {
        for (auto vertex : ORIGINAL_VERTICES) {
            if (sqrt(vertex[0] * vertex[0] +
                     vertex[1] * vertex[1] +
                     vertex[2] * vertex[2]) != 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * indexes of the coordinates of the 20 triangles
     * order of vertices is right-hand-rule
     * counter-clockwise looking at sphere from outside
     */
    const TriangleList ORIGINAL_TRIANGLES = {
            {0, 2, 1},
            {0, 3, 2},
            {0, 4, 3},
            {0, 5, 4},
            {0, 1, 5},
            {1, 2, 6},
            {2, 3, 7},
            {3, 4, 8},
            {4, 5, 9},
            {5, 1, 10},
            {6, 2, 7},
            {8, 4, 9},
            {9, 5, 10},
            {10, 1, 6},
            {7, 3, 8},
            {6, 7, 11},
            {7, 8, 11},
            {8, 9, 11},
            {9, 10, 11},
            {10, 6, 11}
    };

    void geodesic_sphere(VertexList& new_vertices, TriangleList& new_triangles, const double& radius = 1, const unsigned int& resolution = 4) {
        new_vertices.clear();
        new_triangles.clear();

        for (auto triangleI = ORIGINAL_TRIANGLES.begin(); triangleI != ORIGINAL_TRIANGLES.end(); ++triangleI) {
            // assuming this orientation
            //   0
            //  / \
            // 1---2
            Vertex dir_of_hor_lines = (ORIGINAL_VERTICES[triangleI->at(2)] -
                                       ORIGINAL_VERTICES[triangleI->at(1)]) / resolution;

            Vertex dir_of_0_1 = (ORIGINAL_VERTICES[triangleI->at(1)] -
                                 ORIGINAL_VERTICES[triangleI->at(0)]) / resolution;

            Vertex starting_point = ORIGINAL_VERTICES[triangleI->at(0)];
            new_vertices.push_back(starting_point);
            for (int i = 1; i <= resolution; ++i) {
                starting_point += dir_of_0_1;
                new_vertices.push_back(starting_point);

                // now go over on the horizontal line
                Vertex moving_hor_point = starting_point;
                for (int j = 1; j <= i; ++j) {
                    moving_hor_point += dir_of_hor_lines;
                    new_vertices.push_back(moving_hor_point);

                    // every time I move over on this horizontal line
                    // I've completed 2 new triangles if not on the last j
                    // 1 new triangle on the last j
                    size_t li = new_vertices.size() - 1;  // last index
                    // the triangle finished every time
                    new_triangles.push_back({li-1, li, li-(1+i)});  // clockwise
                    // the triangle finished if not on the last j
                    if (j != i) {
                        new_triangles.push_back({li, li-i, li-(i+1)});  // clockwise
                    }
                }  // done moving over on horizontal line
            }  // done with all horizontal lines
        }  // done with all triangles

        // push out all new vertices to radius from origin
        for (auto vertexI = new_vertices.begin(); vertexI != new_vertices.end(); ++vertexI) {
            double mag_of_this_vector = sqrt(vertexI->at(0) * vertexI->at(0) +
                                             vertexI->at(1) * vertexI->at(1) +
                                             vertexI->at(2) * vertexI->at(2));
            *vertexI = *vertexI * (radius / mag_of_this_vector);
        }

        // make sure close vertices are equal
        // threshold is 1/3 of distance between first two vertices
        double threshold = distance_between_2_vertices(new_vertices[0], new_vertices[1]) / 3.0;
        for (auto vertexI = new_vertices.begin(); vertexI != new_vertices.end(); ++vertexI) {
            auto vertexJ = vertexI;
            ++vertexJ;
            while (vertexJ != new_vertices.end()) {
                if (distance_between_2_vertices(*vertexI, *vertexJ) < threshold) {
                    *vertexJ = *vertexI;
                }
                ++vertexJ;
            }
        }

        // TODO: consolidate duplicate vertices
    }

    void create_stl(const std::string& filename, const VertexList& vertices, const TriangleList& triangles) {
        std::ofstream file_out(filename);

        if (! file_out)
            return;

        file_out << "solid geodesic_sphere\n";

        for (auto triangleI = triangles.begin(); triangleI != triangles.end(); ++triangleI) {
            file_out << "  facet normal 0 0 0\n    outer loop\n";
            for (auto indexI = triangleI->begin(); indexI != triangleI->end(); ++indexI) {
                file_out << "      vertex " << vertices[*indexI][0] << ' ' <<
                                               vertices[*indexI][1] << ' ' <<
                                               vertices[*indexI][2] << '\n';
            }
            file_out << "    endloop\n  endfacet\n";
        }

        file_out << "endsolid geodesic_sphere\n";
    }
};

#endif //GEODESICSPHERE_GS_H
