#ifndef BINOBJECT_H
#define BINOBJECT_H

#include <string>
#include <fstream>
#include <vector>
#include <chrono>

// DGtal
#include <DGtal/base/Common.h>                      // We must always include this for using the DGtal: Z3i
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/images/ImageContainerBySTLVector.h> // For setFromFile2D
#include <DGtal/images/ImageSelector.h>             // For setFromFile2D
// PHAT
#include "../include/phat/compute_persistence_pairs.h" // wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm

using namespace DGtal;

struct TB_pair_info
{
    int dimension = 0;
    double thickness = 0;
    std::vector<int> t_center;
    double breadth = 0;
    std::vector<int> b_center;
};

/**
 * @brief The BinObject class
 * Class for the binary object (2D image or 3D volume) and the computation of
 * its thickness-breadth pairs
 */
class BinObject
{
public:
    typedef ImageSelector<Z2i::Domain, int>::Type Image2D;
    typedef ImageSelector<Z3i::Domain, int>::Type Image3D;

    BinObject(const std::string &filename, int dimension);
    void set_benchmark();

    void set_from_pgm();
    void set_from_pgm_2D();
    void set_from_pgm_3D();

    void compute_persistence();
    void write_tb_pairs(const std::vector<TB_pair_info> &tb_pairs) const;
    std::vector<int> coordinates_voxel(int index) const;
    std::vector<int> coordinates_cube(int index) const;
    int index_of_cube(int x, int y, int z) const;
    int index_of_voxel(int x, int y, int z) const;
    int dim_of_cube(int cube) const;
    std::list<int> all_faces(int voxel) const;
    void faces(int cube, std::list<int> &f) const;
    std::list<int> incident_voxels(int cube) const;
    double value_of_cube(int cube) const;
    int voxel_of_cube(int cube) const;

    double elapsed_sec(std::chrono::steady_clock::time_point start_time) const;

private:
    std::string m_filename;         /// name of the file containing the object
    int m_dimension;                /// dimension of the object (2 or 3)
    std::vector<unsigned> m_size;   /// number of voxels in the object along each axis
    std::vector<double> m_sdt;      /// signed distance transform
    bool m_benchmark;               /// output running times in one line
};

#endif // BINOBJECT_H
