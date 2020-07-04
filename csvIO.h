#pragma once

void writeToCsv(std::string filename, std::vector<std::pair<std::string, std::vector<double>>> dataset);
std::vector<std::pair<std::string, std::vector<double>>> readFromCsv(std::string filename);