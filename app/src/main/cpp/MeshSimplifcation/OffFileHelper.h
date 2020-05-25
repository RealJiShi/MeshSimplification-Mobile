#pragma once

#include <fstream>
#include <sstream>
#include <vector>

class OffFileHelper
{
public:
    static bool load(const std::string &file_path,
                     std::vector<double> &vertices,
                     std::vector<uint16_t> &indices)
    {
        std::ifstream fin(file_path);
        if (!fin)
        {
            return false;
        }

        // get "OFF"
        std::string line;
        while (std::getline(fin, line))
        {
            if (line.length() > 1 && line[0] != '#')
            {
                if (line.find("OFF") != std::string::npos)
                {
                    break;
                }
            }
        }

        // get vertex number and face number
        unsigned int nVert = 0;
        unsigned int nFace = 0;
        unsigned int nEdge = 0;
        while (std::getline(fin, line))
        {
            if (line.length() > 1 && line[0] != '#')
            {
                std::istringstream info;
                info.str(line);
                info >> nVert >> nFace >> nEdge;
                break;
            }
        }

        // get vertex
        vertices.resize(nVert * 3);
        for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert)
        {
            std::getline(fin, line);
            std::stringstream info;
            info.str(line);
            info >> vertices[i_vert * 3] >> vertices[i_vert * 3 + 1] >> vertices[i_vert * 3 + 2];
        }

        // get face
        indices.resize(nFace * 3);
        for (unsigned int i_face = 0; i_face < nFace; ++i_face)
        {
            std::getline(fin, line);
            std::stringstream info;
            info.str(line);
            uint16_t _temp;
            info >> _temp >> indices[i_face * 3] >> indices[i_face * 3 + 1] >> indices[i_face * 3 + 2];
        }
        fin.close();

        //
        return true;
    }

    static bool store(const std::string &file_path,
                      const std::vector<double> &vertices,
                      const std::vector<uint16_t> &indices)
    {
        std::ofstream fout(file_path);
        if (!fout)
        {
            return false;
        }

        unsigned int nVert = static_cast<unsigned int>(vertices.size() / 3);
        unsigned int nFace = static_cast<unsigned int>(indices.size() / 3);

        fout << "OFF\n";
        fout << nVert << " " << nFace << " " << 0 << std::endl;
        for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert)
        {
            fout << vertices[i_vert * 3] << " " << vertices[i_vert * 3 + 1] << " " << vertices[i_vert * 3 + 2] << std::endl;
        }

        for (unsigned int i_face = 0; i_face < nFace; ++i_face)
        {
            fout << 3 << " " << indices[i_face * 3] << " " << indices[i_face * 3 + 1] << " " << indices[i_face * 3 + 2] << std::endl;
        }
        fout.close();

        return true;
    }
};