
#include <iostream>
#include <FileHandler.h>
#include <vector>
#include <string>
#include <filesystem>

std::vector<std::string> FileHandler::getDirectoryFiles()
{
    // List of file names.
    std::vector<std::string> files{};
    std::string file_prefix = "/home/antu/Python_Projects/graph_research/";
    const fs::path target_path{file_prefix};

    try {
        for (auto const& dir_entry : fs::directory_iterator{target_path})
        {
            if (fs::is_regular_file(dir_entry.path()))
            {
                files.push_back(file_prefix + dir_entry.path().filename().string());
            }
        }
    }
    catch (fs::filesystem_error const& ex)
    {
        std::cout << "Error occurred during file operations!\n" << ex.what() << std::endl;
        return files;
    }
    return files;
}