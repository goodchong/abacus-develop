#ifndef GLOBAL_FILE_H
#define GLOBAL_FILE_H

#include <string>

namespace ModuleBase
{
namespace Global_File
{

void make_h0_output_dir(int rank,
                        const std::string& output_dir,
                        const std::string& log_file);
void close_all_log(int rank);

} // namespace Global_File
} // namespace ModuleBase

#endif
