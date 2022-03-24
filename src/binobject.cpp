#include "binobject.h"

// DGtal
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/io/boards/Board2D.h"
#include <DGtal/shapes/Shapes.h>

// PHAT
#include "../include/phat/compute_persistence_pairs.h" // wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include "../include/phat/representations/default_representations.h" // main data structure (choice affects performance)
#include "../include/phat/algorithms/standard_reduction.h" // algorithm (choice affects performance)
#include "../include/phat/algorithms/chunk_reduction.h"
#include "../include/phat/algorithms/row_reduction.h"
#include "../include/phat/algorithms/twist_reduction.h"

#if defined(__APPLE__)
// Boost
#include <boost/algorithm/minmax_element.hpp>
#endif

#include <string>

BinObject::BinObject(const std::string &filename, int dimension)
    : m_filename(filename),
      m_dimension(dimension),
      m_benchmark(false)
{
}


void BinObject::set_benchmark()
{
    m_benchmark = true;
    std::cout << "Filename & voxels & sdt & matrix & pairs" << std::endl;
}


/**
 * @brief Read a PGM file
 */
void BinObject::set_from_pgm()
{
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    if (m_dimension == 2)
        set_from_pgm_2D();
    else
        set_from_pgm_3D();

    if (m_benchmark)
    {
        std::cout << m_filename << " & "
                  << m_size[0]*m_size[1]*m_size[2] << " & "
                  << elapsed_sec(start_time) << " & ";
    }
    else
        std::clog << "Signed distance transform computed: " << elapsed_sec(start_time) << std::endl;
}

/**
 * @brief Reads the PGM file and compute the signed distance transform
 */
void BinObject::set_from_pgm_2D()
{
    assert(m_dimension == 2);

    const Image2D image = DGtal::GenericReader<Image2D>::import(m_filename + ".pgm"); // import the image
    const Z2i::Domain dom = image.domain();   // define its domain for knowing the bounds

    std::cout << "Name: " << m_filename + ".pgm" << std::endl;
    std::cout << "Size of the image: "
              << dom.myUpperBound[0] + 1 << " x "
              << dom.myUpperBound[1] + 1 << std::endl;

    m_size.resize(m_dimension);
    for (int i = 0; i < m_dimension; i++)
        m_size[i] = dom.myUpperBound[i] + 1;

//    for (Z2i::Domain::Iterator it = dom.begin(); it != dom.end(); it++)   // printing the values of the image
//        std::cout << *it << " value = " << image(*it) << std::endl;
//    std::cout << "Image loaded from <" << m_filename << ">." << std::endl;

    typedef functors::IntervalForegroundPredicate<Image2D> Binarizer;   // for binarizing
    Binarizer b_in(image, 0, 255);    // the interval ]0,255] is considered true
    typedef DGtal::DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;   // the Distance Transform for L2 (Euclidean)
    DTL2 dt_in(&image.domain(), &b_in, &Z2i::l2Metric);   // compute the distance transform in dt_in

    Binarizer b_out(image, -1, 0);    // b_out are the 0-valued pixels
    DTL2 dt_out(&image.domain(), &b_out, &Z2i::l2Metric);

    m_sdt.resize(dom.size());
    std::size_t i = 0;
    for (Z2i::Domain::ConstIterator it = dom.begin(); it != dom.end(); it++)
    {
        m_sdt.at(i) = dt_out(*it) - dt_in(*it);
        i++;
    }
}

/**
 * @brief Reads the PGM (3D) file
 * Compute the distance transform of the foreground and the background
 * Write the shifted (so the values start by 0) SDT of the picture
 */
void BinObject::set_from_pgm_3D()
{
    assert(m_dimension == 3);

    const Image3D volume = DGtal::GenericReader<Image3D>::import(m_filename + ".pgm");
    const Z3i::Domain dom = volume.domain();

    if (!m_benchmark)
    {
        std::cout << "Name: " << m_filename + ".pgm" << std::endl;
        std::cout << "Size of the image: "
                  << dom.myUpperBound[0] + 1 << " x "
                  << dom.myUpperBound[1] + 1 << " x "
                  << dom.myUpperBound[2] + 1
                  << " [" << dom.size() << "]" << std::endl;
    }

    m_size.resize(m_dimension);
    for (int i = 0; i < m_dimension; i++)
        m_size[i] = dom.myUpperBound[i] + 1;

//    for (Z3i::Domain::Iterator it = dom.begin(); it != dom.end(); it++) // print the values of the voxels for debugging
//        std::cout << *it << " value = " << volume(*it) << std::endl;
//    std::cout << "Volume loaded from " << fileName << "." << std::endl;

    typedef functors::IntervalForegroundPredicate<Image3D> Binarizer;
    const Binarizer b_in(volume, 0, 255);    // the interval ]0,255] is considered true
    typedef DGtal::DistanceTransformation<Z3i::Space, Binarizer, Z3i::L2Metric> DTL2;   // the Distance Transform for L2 (Euclidean)
    const DTL2 dt_in(&volume.domain(), &b_in, &Z3i::l2Metric);   // compute the distance transform in dt_in

    const Binarizer b_out(volume, -1, 0);
    const DTL2 dt_out(&volume.domain(), &b_out, &Z3i::l2Metric);

    m_sdt.resize(dom.size());
    std::size_t i = 0;
    for (Z3i::Domain::ConstIterator it = dom.begin(); it != dom.end(); it++)
    {
        m_sdt.at(i) = dt_out(*it) - dt_in(*it);
        i++;
    }
}

/**
 * @brief BinObject::compute_persistence Computes the persistent homology with PHAT
 * @param benchmark if true, output the running times in one line
 */
void BinObject::compute_persistence()
{
    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    // sort voxels by signed distance transform value
    std::vector<std::pair<double, int>> voxels(m_sdt.size());
    for (std::size_t i = 0; i < voxels.size(); i++)
        voxels.at(i) = std::make_pair(m_sdt.at(i), i);
    std::sort(voxels.begin(), voxels.end());
    if (!m_benchmark)
        std::clog << "voxels sorted: " << elapsed_sec(start_time) << std::endl;

    // build filtered cubical complex
    const std::size_t nb_cubes = (2*m_size[0]+1)*(2*m_size[1]+1)*(2*m_size[2]+1);
    std::vector<int> pos_of_cube(nb_cubes, -1); // position of a cube (given its index) in the filter
    std::vector<int> cube_at_pos(nb_cubes, -1); // the inverse: the (index of the) cube at a given position in the filter
    std::vector<double> sdt_cube(nb_cubes); // the function value on the cubes
    int filter = 0;
    for (auto it = voxels.cbegin(); it != voxels.cend(); ++it)
    {
        // add faces of this voxels, sorted by dimension, unless they are already added
        const std::list<int> cubes = all_faces(it->second);
        for (int cube : cubes)
        {
            if (pos_of_cube.at(cube) < 0) // if that cube is not already added
            {
                pos_of_cube.at(cube) = filter;
                cube_at_pos.at(filter) = cube;
                sdt_cube.at(cube) = it->first;
                filter++;
            }
        }
    }
    if (!m_benchmark)
        std::clog << "filter built: " << elapsed_sec(start_time) << std::endl;

    // make the boundary matrix
    // first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix<phat::bit_tree_pivot_column> boundary_matrix;
    boundary_matrix.set_num_cols(nb_cubes);
    std::vector<phat::index> temp_col; // current column of the matrix
    for (std::size_t i = 0; i < nb_cubes; i++)
    {
        boundary_matrix.set_dim(i, dim_of_cube(cube_at_pos.at(i)));
        std::list<int> f; // list of faces of the cube
        faces(cube_at_pos.at(i), f);
        temp_col.resize(f.size());
        int j = 0;
        for (int face : f)
        {
            temp_col.at(j) = pos_of_cube.at(face);
            j++;
        }
        std::sort(temp_col.begin(), temp_col.end()); // @note This is very important!
        boundary_matrix.set_col(i, temp_col);
        f.clear();
    }
    if (m_benchmark)
        std::cout << elapsed_sec(start_time) << " & ";
    else
    {
        std::clog << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns and " << boundary_matrix.get_num_entries() << " entries." << std::endl;
        std::clog << "Boundary matrix built: " << elapsed_sec(start_time) << std::endl;
    }
    start_time = std::chrono::steady_clock::now();
//    boundary_matrix.save_ascii("phat_boundary.dat"); // write the boundary matrix for debugging

    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs<phat::chunk_reduction>(pairs, boundary_matrix);
    if (!m_benchmark)
        std::clog << "Persistent homology computed: " << elapsed_sec(start_time) << std::endl;
//    pairs.sort();
//    pairs.save_ascii("pairs.dat");
//    std::cout << "There are " << pairs.get_num_pairs() << " persistence pairs: " << std::endl;

    // get the thickness-breadth pairs
    std::vector<TB_pair_info> tb_pairs;
    for (phat::index idx = 0; idx < pairs.get_num_pairs(); idx++)
    {
//        std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
        const int t_cell = cube_at_pos.at(pairs.get_pair(idx).first);
        const int b_cell = cube_at_pos.at(pairs.get_pair(idx).second);
        const double t_value = sdt_cube.at(t_cell);
        const double b_value = sdt_cube.at(b_cell);
        if (t_value < 0 && b_value > 0)
        {
            const std::vector<int> t_ball = coordinates_voxel(voxel_of_cube(t_cell));
            const std::vector<int> b_ball = coordinates_voxel(voxel_of_cube(b_cell));
//            const std::vector<int> t_ball = coordinates_cube(t_cell); // @warning I am getting just any voxel adjacent to this cell!
//            const std::vector<int> b_ball = coordinates_cube(b_cell);
            TB_pair_info cur_tb_pair;
            cur_tb_pair.dimension = dim_of_cube(t_cell);
            cur_tb_pair.thickness = -t_value;
            cur_tb_pair.t_center = t_ball;
//            for (auto x : t_ball)
//                cur_tb_pair.t_center.push_back((x-1)/2);
            cur_tb_pair.breadth = b_value;
            cur_tb_pair.b_center = b_ball;
//            for (auto x : b_ball)
//                cur_tb_pair.b_center.push_back((x-1)/2);
            tb_pairs.push_back(cur_tb_pair);
        }
    }
    TB_pair_info last_tb_pair;
    last_tb_pair.dimension = 0;
    last_tb_pair.thickness = -voxels.front().first;
    for (auto x : coordinates_cube(cube_at_pos.front()))
        last_tb_pair.t_center.push_back((x-1)/2);
    last_tb_pair.breadth = -1;
    tb_pairs.push_back(last_tb_pair);

    write_tb_pairs(tb_pairs);

    if (m_benchmark)
        std::cout << elapsed_sec(start_time) << " \\\\" << std::endl;
    else
        std::clog << "TB-pairs computed: " << elapsed_sec(start_time) << std::endl;
}

void BinObject::write_tb_pairs(const std::vector<TB_pair_info> &tb_pairs) const
{
    std::ofstream file(m_filename + ".json", std::ios::out | std::ios::trunc);
    if (!(file))
    {
        std::cerr << "Error in BinObject::write_tb_pairs(): impossible to create the output text file." << std::endl;
        exit(EXIT_FAILURE);
    }

    file << "{" << std::endl;
    file << "\t\"meta\": {" << std::endl;
    file << "\t\t\"filename\": \"" << m_filename << ".pgm\"" << std::endl;
    file << "\t}," << std::endl;
    file << "\t\"pairs\": [" << std::endl;

    for (auto it = tb_pairs.cbegin(); it != tb_pairs.cend(); ++it) {
        file << "\t\t{" << std::endl;
        file << "\t\t\t\"dimension\": " << it->dimension << "," << std::endl;
        file << "\t\t\t\"thickness\": " << it->thickness << "," << std::endl;
        file << "\t\t\t\"thickness_center\": [" << it->t_center[0] << ", " << it->t_center[1] << ", " << it->t_center[2] << "], " << std::endl;

        if (it->breadth >= 0) {
            file << "\t\t\t\"breadth\": " << it->breadth << "," << std::endl;
            file << "\t\t\t\"breadth_center\": [" << it->b_center[0] << ", " << it->b_center[1] << ", " << it->b_center[2] << "]" << std::endl;
            file << "\t\t}," << std::endl;
        }
        else
        {
            file << "\t\t\t\"breadth\": " << it->breadth << std::endl;
            file << "\t\t}" << std::endl;
        }
    }
    file << "\t]" << std::endl;
    file << "}" << std::endl;
    file.close();
}

/**
 * @brief BinObject::index_of_cube
 * @return The index of a cube in a cubical complex given its coordinates
 */
int BinObject::index_of_cube(int x, int y, int z) const
{
    assert(0 <= x && x < (2*m_size[0]+1));
    assert(0 <= y && y < (2*m_size[1]+1));
    assert(0 <= z && z < (2*m_size[2]+1));

    int index = 0;
    index += x;
    index += y * (2*m_size[0]+1);
    index += z * (2*m_size[0]+1)*(2*m_size[1]+1);
    return index;
}

/**
 * @brief BinObject::index_of_voxel
 * @return The index of voxel (x,y,z) or -1 if it is outside the box
 */
int BinObject::index_of_voxel(int x, int y, int z) const
{
    if (x < 0 || x >= (int)m_size.at(0)) return -1;
    if (z < 0 || y >= (int)m_size.at(1)) return -1;
    if (x < 0 || z >= (int)m_size.at(2)) return -1;

    int index = 0;
    index += x;
    index += y * m_size[0];
    index += z * m_size[0]*m_size[1];
    return index;
}

/**
 * @brief BinObject::dim_of_cube
 * @param cube
 * @return The dimension of the cube given its index
 */
int BinObject::dim_of_cube(int cube) const
{
    const std::vector<int> coord = coordinates_cube(cube);
    return coord[0]%2 + coord[1]%2 + coord[2]%2;
}

/**
 * @brief BinObject::faces
 * @return The indices of the faces of the 3-cube associated to the voxels with
 * index @param voxel. The cubes are sorted by dimension
 *
 * @todo Do not use the function index_of_cube and compute the indices locally
 */
std::list<int> BinObject::all_faces(int voxel) const
{
    const std::vector<int> coord = coordinates_voxel(voxel); // voxel coordinates
    const int x = coord[0];
    const int y = coord[1];
    const int z = coord[2];

    std::list<int> f;
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z  )); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z  )); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z  )); // dim 2
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y  , 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y  , 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y  , 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x  , 2*y+2, 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x  , 2*y+1, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z+2)); // dim 0
    f.push_back(index_of_cube(2*x+2, 2*y+2, 2*z+1)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+2, 2*y+1, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z+2)); // dim 1
    f.push_back(index_of_cube(2*x+1, 2*y+2, 2*z+1)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z+2)); // dim 2
    f.push_back(index_of_cube(2*x+1, 2*y+1, 2*z+1)); // dim 3
    return f;
}

/**
 * @brief BinObject::coordinates_voxel
 * @param index
 * @return The coordinates of a voxel in a discrete object given its index
 */
std::vector<int> BinObject::coordinates_voxel(int index) const
{
    std::vector<int> coord(3);
    coord[0] = index % m_size[0];
    index /= m_size[0];
    coord[1] = index % m_size[1];
    index /= m_size[1];
    coord[2] = index;
    return coord;
}

/**
 * @brief BinObject::coordinates_cube
 * @param index
 * @return The coordinates of a cube in a cubical complex given its index
 */
std::vector<int> BinObject::coordinates_cube(int index) const
{
    std::vector<int> coord(3);
    coord[0] = index % (2*m_size[0]+1);
    index /= (2*m_size[0]+1);
    coord[1] = index % (2*m_size[1]+1);
    index /= (2*m_size[1]+1);
    coord[2] = index;
    return coord;
}

/**
 * @brief BinObject::faces
 * @param cube_index
 * @return The list of the faces of @param cube_index
 */
void BinObject::faces(int cube, std::list<int> &f) const
{
    const std::vector<int> coord = coordinates_cube(cube);
    int offset = 1;
    for (int i = 0; i < 3; i++)
    {
        if (coord.at(i) % 2 == 1)
        {
            f.push_back(cube - offset);
            f.push_back(cube + offset);
        }
        offset *= (2*m_size[i]+1);
    }
}


/**
 * @brief BinObject::incident_voxels
 * @param cube
 * @return The list of the indices of the voxels that are adjacent to a given cube.
 * @note I do not use this anymore, but I keep it because it was very long to write.
 */
std::list<int> BinObject::incident_voxels(int cube) const
{
    std::list<int> voxels;
    const std::vector<int> coord = coordinates_cube(cube);
    const int x = coord[0];
    const int y = coord[1];
    const int z = coord[2];
    if (x%2 == 0 && y%2 == 0 && z%2 == 0)
    {
        voxels.push_back(index_of_voxel((x+1)/2, (y+1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y+1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x+1)/2, (y-1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y-1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x+1)/2, (y+1)/2, (z-1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y+1)/2, (z-1)/2));
        voxels.push_back(index_of_voxel((x+1)/2, (y-1)/2, (z-1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y-1)/2, (z-1)/2));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 0)
    {
        voxels.push_back(index_of_voxel((x  )/2, (y+1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x  )/2, (y-1)/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x  )/2, (y+1)/2, (z-1)/2));
        voxels.push_back(index_of_voxel((x  )/2, (y-1)/2, (z-1)/2));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 0)
    {
        voxels.push_back(index_of_voxel((x+1)/2, (y  )/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y  )/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x+1)/2, (y  )/2, (z-1)/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y  )/2, (z-1)/2));
    }
    else if (x%2 == 0 && y%2 == 0 && z%2 == 1)
    {
        voxels.push_back(index_of_voxel((x+1)/2, (y+1)/2, (z  )/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y+1)/2, (z  )/2));
        voxels.push_back(index_of_voxel((x+1)/2, (y-1)/2, (z  )/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y-1)/2, (z  )/2));
    }
    else if (x%2 == 1 && y%2 == 1 && z%2 == 0)
    {
        voxels.push_back(index_of_voxel((x  )/2, (y  )/2, (z+1)/2));
        voxels.push_back(index_of_voxel((x  )/2, (y  )/2, (z-1)/2));
    }
    else if (x%2 == 1 && y%2 == 0 && z%2 == 1)
    {
        voxels.push_back(index_of_voxel((x  )/2, (y+1)/2, (z  )/2));
        voxels.push_back(index_of_voxel((x  )/2, (y-1)/2, (z  )/2));
    }
    else if (x%2 == 0 && y%2 == 1 && z%2 == 1)
    {
        voxels.push_back(index_of_voxel((x+1)/2, (y  )/2, (z  )/2));
        voxels.push_back(index_of_voxel((x-1)/2, (y  )/2, (z  )/2));
    }
    else
    {
        voxels.push_back(index_of_voxel((x  )/2, (y  )/2, (z  )/2));
    }
    return voxels;
}





/**
 * @brief BinObject::value_of_cube
 * @param cube
 * @return Given the index of a cube, get its incident voxels (3-cubes) and get the minimum value
 * @note This is not used anymore
 */
double BinObject::value_of_cube(int cube) const
{
    const std::list<int> voxels = incident_voxels(cube);
    double min_sdt = std::numeric_limits<double>::max();
    for (int voxel : voxels)
        if (voxel >= 0)
            min_sdt = std::min(min_sdt, m_sdt.at(voxel));
    assert(min_sdt != std::numeric_limits<double>::max());
    return min_sdt;
}


/**
 * @brief BinObject::voxel_of_cube
 * @param cube
 * @return Given a cube in the filter, return the index of the voxels that adds
 * it to the filter. That is, the adjacent voxel with least signed distance
 * transform value
 */
int BinObject::voxel_of_cube(int cube) const
{
    const std::list<int> voxels = incident_voxels(cube);
    double min_sdt = std::numeric_limits<double>::max();
    int first_voxel = 0;
    for (int voxel : voxels)
        if (voxel >= 0 && m_sdt.at(voxel) < min_sdt)
        {
            min_sdt = m_sdt.at(voxel);
            first_voxel = voxel;
        }
    assert(min_sdt != std::numeric_limits<double>::max());
    return first_voxel;
}


double BinObject::elapsed_sec(std::chrono::steady_clock::time_point start_time) const
{
    auto cur_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = cur_time - start_time;
    return elapsed_seconds.count();
}
