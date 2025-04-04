#ifndef SUBGRAPHMATCHINGMAIN_FILEHANDLER_H
#define SUBGRAPHMATCHINGMAIN_FILEHANDLER_H

#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;


class FileHandler {

public:
    static std::vector<std::string> getDirectoryFiles();
    

};


#endif //SUBGRAPHMATCHINGMAIN_FILEHANDLER_H