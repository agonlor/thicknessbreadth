#include "binobject.h"

#include <iostream>
#include <sstream>  // convert argv to int
#include <string>

void help(int argc, char **argv)
{
    std::cout << " == Thickness-breadth ==" << std::endl;
    std::cout << " Compute the thickness-breadth pairs of a discrete object" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " <input file name> [OPTIONS]" << std::endl;
    std::cout << "  <input file name>: a pgm file of type P2 (ASCII)" << std::endl;
    std::cout << "OPTIONS" << std::endl;
    std::cout << "    -d: dimension of the pgm file" << std::endl;
    std::cout << "        it can be 2 or 3 (by default)" << std::endl;
    std::cout << "    --benchmark: output running times" << std::endl;
}

void usage(int argc, char **argv)
{
    std::cout << "Usage: " << argv[0] << " <input file name> [OPTIONS]" << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "Example: " << argv[0] << " volume.pgm" << std::endl;
    std::cout << "More details: " << argv[0] << " -h" << std::endl;
}


int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "-h"))
    {
        help(argc, argv);
        return 0;
    }
    if (argc <  2)
    {
        usage(argc, argv);
        return 0;
    }

    std::string fileName = argv[1];
    fileName.erase(fileName.end()-4, fileName.end()); // filename is argvc[1] without ".pgm"

    int dimension = 3;
    bool benchmark = false;
    for (int i = 2; i < argc; i++)
    {
        std::string arg = argv[i];
        /* dimension */
        if (arg.compare("-d") == 0)
        {
            if (i+1 < argc) {
                arg = argv[++i];
                std::istringstream iss(argv[i]); iss >> dimension;
                if (dimension < 2 || dimension > 3)
                {
                    std::cerr << "Error: dimension (-d) must be 2 or 3." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else { std::cerr << argv[i] << " -d option needs an int value." << std::endl; exit(EXIT_FAILURE); }
        }
        else if (arg.compare("--benchmark") == 0)
        {
            benchmark = true;
        }
    }

    BinObject I(fileName, dimension);
    if (benchmark)
        I.set_benchmark();

    I.set_from_pgm();
    I.compute_persistence();

    return EXIT_SUCCESS;
}
