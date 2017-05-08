// library to create geodesic sphere

#ifndef GEODESICSPHERE_GS_H
#define GEODESICSPHERE_GS_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>

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

    void geodesic_sphere(VertexList& output_vertices,
                         TriangleList& output_triangles,
                         const double& radius = 1,
                         const unsigned int& resolution = 4)
    {
        output_vertices.clear();
        output_triangles.clear();

        for (auto triangleI = ORIGINAL_TRIANGLES.begin(); triangleI != ORIGINAL_TRIANGLES.end(); ++triangleI) {
            // assuming this orientation
            //   0
            //  / \
            // 1---2

            // direction of horizontal lines
            Vertex dir_of_hor_lines = (ORIGINAL_VERTICES[triangleI->at(2)] -
                                       ORIGINAL_VERTICES[triangleI->at(1)]) / resolution;

            // direction of the line from 0 to 1
            Vertex dir_of_0_1 = (ORIGINAL_VERTICES[triangleI->at(1)] -
                                 ORIGINAL_VERTICES[triangleI->at(0)]) / resolution;

            Vertex starting_point = ORIGINAL_VERTICES[triangleI->at(0)];
            output_vertices.push_back(starting_point);

            // for each horizontal line
            for (int i = 1; i <= resolution; ++i) {
                starting_point += dir_of_0_1;
                output_vertices.push_back(starting_point);

                // now move right along the horizontal line
                Vertex moving_hor_point = starting_point;
                for (int j = 1; j <= i; ++j) {
                    moving_hor_point += dir_of_hor_lines;
                    output_vertices.push_back(moving_hor_point);

                    // every time I move over on this horizontal line
                    // I've completed 2 new triangles if not on the last j
                    // 1 new triangle on the last j
                    size_t li = output_vertices.size() - 1;  // last index
                    // the triangle finished every time
                    output_triangles.push_back({li-1, li, li-(1+i)});  // counter-clockwise
                    // the triangle finished if not on the last j
                    if (j != i) {
                        output_triangles.push_back({li, li-i, li-(i+1)});  // counter-clockwise
                    }
                }  // done moving over on horizontal line
            }  // done with all horizontal lines
        }  // done with all triangles

        // push out all new vertices to radius from origin
        for (auto vertexI = output_vertices.begin(); vertexI != output_vertices.end(); ++vertexI) {
            double mag_of_this_vector = sqrt(vertexI->at(0) * vertexI->at(0) +
                                             vertexI->at(1) * vertexI->at(1) +
                                             vertexI->at(2) * vertexI->at(2));
            *vertexI = *vertexI * (radius / mag_of_this_vector);
        }

        // to make sure there are no cracks between triangles
        // make sure close vertices are equal
        // (floating point arithmetic might have made them unequal)
        // threshold is 1/3 of distance between first two vertices
        double threshold = distance_between_2_vertices(output_vertices[0], output_vertices[1]) / 3.0;
        for (auto vertexI = output_vertices.begin(); vertexI != output_vertices.end(); ++vertexI) {
            auto vertexJ = vertexI;
            ++vertexJ;
            while (vertexJ != output_vertices.end()) {
                if (distance_between_2_vertices(*vertexI, *vertexJ) < threshold) {
                    *vertexJ = *vertexI;
                }
                ++vertexJ;
            }
        }

        // TODO: consolidate duplicate vertices
    }

    Vertex calculate_facet_normal(const Triangle& triangle, const VertexList& vertices) {
        // https://math.stackexchange.com/questions/305642/how-to-find-surface-normal-of-a-triangle/
        Vertex v = vertices[triangle[1]] - vertices[triangle[0]];
        Vertex w = vertices[triangle[2]] - vertices[triangle[0]];
        return {
            (v[1] * w[2]) - (v[2] * w[1]),
            (v[2] * w[0]) - (v[0] * w[2]),
            (v[0] * w[1]) - (v[1] * w[0])
        };
    }

    double get_translator(bool allow_negative_coordinates, const VertexList& vertices) {
        return (int)(! allow_negative_coordinates) * distance_between_2_vertices(vertices[0], {0, 0, 0});
    }

    bool ends_with_case_insensitive(const std::string& original, const std::string& ending)
    {
        if (ending.size() > original.size()) return false;
        auto original_ending_it = original.rbegin();
        for (auto ending_it = ending.rbegin(); ending_it != ending.rend(); ++ending_it) {
            if (tolower(*original_ending_it) != tolower(*ending_it)) {
                return false;
            }
            ++original_ending_it;
        }
        return true;
    }

    std::string confirm_valid_filename(const std::string& filename, const std::string& extension_without_dot) {
        if (filename.empty()) {
            return "untitled." + extension_without_dot;
        }
        if (filename.size() < extension_without_dot.size() + 2 ||
            filename[filename.size() - (extension_without_dot.size() + 1)] != '.' ||
            (! ends_with_case_insensitive(filename, extension_without_dot)))
        {
            return filename + '.' + extension_without_dot;
        }
        return filename;
    }

    void create_ascii_stl(std::string filename,
                          const VertexList& vertices,
                          const TriangleList& triangles,
                          bool allow_negative_coordinates = false)
    {
        // verify filename ends with .stl
        filename = confirm_valid_filename(filename, "stl");

        // in case we need to translate coordinates
        double translator = get_translator(allow_negative_coordinates, vertices);

        std::ofstream file_out(filename);

        if (! file_out)
            exit(1);

        file_out << "solid geodesic_sphere\n";

        for (auto triangleI = triangles.begin(); triangleI != triangles.end(); ++triangleI) {

            Vertex n = calculate_facet_normal(*triangleI, vertices);

            file_out << "  facet normal " << n[0] << ' ' << n[1] << ' ' << n[2] << "\n    outer loop\n";
            for (auto indexI = triangleI->begin(); indexI != triangleI->end(); ++indexI) {
                file_out << "      vertex " << vertices[*indexI][0] + translator << ' ' <<
                                               vertices[*indexI][1] + translator << ' ' <<
                                               vertices[*indexI][2] + translator << '\n';
            }
            file_out << "    endloop\n  endfacet\n";
        }

        file_out << "endsolid geodesic_sphere\n";
    }

    // http://stackoverflow.com/questions/33134811/writing-binary-stl-files-in-c
    void create_binary_stl(std::string filename,
                           const VertexList& vertices,
                           const TriangleList& triangles,
                           bool allow_negative_coordinates = false)
    {
        // verify filename ends with .stl
        filename = confirm_valid_filename(filename, "stl");

        // in case we need to translate coordinates
        float translator = (float)get_translator(allow_negative_coordinates, vertices);

        char head[80];
        std::strncpy(head, filename.c_str(), sizeof(head) - 1);
        char attribute[2] = {0, 0};
        uint32_t triangle_count = (uint32_t) triangles.size();

        std::ofstream file_out;

        file_out.open(filename, std::ios::out | std::ios::binary);

        if (!file_out)
            exit(1);

        file_out.write(head, sizeof(head));
        file_out.write((char *) &triangle_count, 4);

        // write down every triangle
        for (auto triangleI = triangles.begin(); triangleI != triangles.end(); ++triangleI) {
            //normal vector coordinates
            Vertex n = calculate_facet_normal(*triangleI, vertices);
            float x = (float) n[0];
            float y = (float) n[1];
            float z = (float) n[2];

            file_out.write((char *) &x, 4);
            file_out.write((char *) &y, 4);
            file_out.write((char *) &z, 4);

            for (auto indexI = triangleI->begin(); indexI != triangleI->end(); ++indexI) {
                x = (float)vertices[*indexI][0] + translator;
                y = (float)vertices[*indexI][1] + translator;
                z = (float)vertices[*indexI][2] + translator;

                file_out.write((char *) &x, 4);
                file_out.write((char *) &y, 4);
                file_out.write((char *) &z, 4);
            }

            file_out.write(attribute, 2);
        }
    }
};

#endif //GEODESICSPHERE_GS_H
