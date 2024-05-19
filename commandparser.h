//
// Created by antu on 5/19/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_COMMANDPARSER_H
#define SUBGRAPHMATCHINGMAIN_COMMANDPARSER_H


#include <string>
#include <algorithm>
#include <vector>
class CommandParser {
private:
    std::vector<std::string> tokens_;

public:
    CommandParser(const int argc, char **argv);
    const std::string getCommandOption(const std::string &option) const;
    bool commandOptionExists(const std::string &option) const;
};

#endif //SUBGRAPHMATCHINGMAIN_COMMANDPARSER_H
