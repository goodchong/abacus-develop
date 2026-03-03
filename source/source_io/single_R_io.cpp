#include "single_R_io.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <string>
#include <charconv>
#include <omp.h>
#include "source_base/timer.h"



static inline void append_sci8(std::string& out, double x)
{
    char buf[64];
    int n = std::snprintf(buf, sizeof(buf), " %.8e", x);
    if (n > 0) out.append(buf, buf + n);
}
static inline void append_sci8(std::string& out, std::complex<double> z)
{
    char buf[160];
    int n = std::snprintf(buf, sizeof(buf), " (%.8e,%.8e)", z.real(), z.imag());
    if (n > 0) out.append(buf, buf + n);
}

inline void write_data(std::ofstream& ofs, const double& data)
{
    ofs << " " << std::fixed << std::scientific << std::setprecision(8) << data;
}
inline void write_data(std::ofstream& ofs, const std::complex<double>& data)
{
    ofs << " (" << std::fixed << std::scientific << std::setprecision(8) << data.real() << ","
        << std::fixed << std::scientific << std::setprecision(8) << data.imag() << ")";
}
template<typename T>
void ModuleIO::output_single_R(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, T>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce)
{
    if (reduce && GlobalV::NPROC > 1)
    {
        output_single_R_reduce(ofs, XR, sparse_threshold, binary, pv);
    }
    else
    {
        if (binary)
            output_single_R_non_reduce_binary(ofs, XR, sparse_threshold, pv);
        else
            output_single_R_non_reduce_txt(ofs, XR, sparse_threshold, pv);
    }
}

template <typename T>
void ModuleIO::output_single_R_reduce(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, T>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv)
{
    T* line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(PARAM.globalv.nlocal + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << PARAM.globalv.global_out_dir << std::to_string(GlobalV::DRANK) + "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    line = new T[PARAM.globalv.nlocal];
    for(int row = 0; row < PARAM.globalv.nlocal; ++row)
    {
        ModuleBase::GlobalFunc::ZEROS(line, PARAM.globalv.nlocal);

        if (pv.global2local_row(row) >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (auto &value : iter->second)
                {
                    line[value.first] = value.second;
                }
            }
        }

        Parallel_Reduce::reduce_all(line, PARAM.globalv.nlocal);

        if (GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < PARAM.globalv.nlocal; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    if (binary)
                    {
                        ofs.write(reinterpret_cast<char*>(&line[col]), sizeof(T));
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        write_data(ofs, line[col]);
                        ofs_tem1 << " " << col;
                    }
                    nonzeros_count++;
                }
            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }
    }

    delete[] line;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), std::ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
            }
        }
        else
        {
            ofs << std::endl;
            ofs_tem1 << std::endl;
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str());
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }
        std::remove(tem1.str().c_str());
    }
}    

template <typename T>
void ModuleIO::output_single_R_non_reduce_txt(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, T>>& XR,
    const double& sparse_threshold,
    const Parallel_Orbitals& pv)
{   
    ModuleBase::timer::tick("ModuleIO", "output_single_R");

    std::vector<char> filebuf(8 * 1024 * 1024 * 8); // 64 MB buffer 
    ofs.rdbuf()->pubsetbuf(filebuf.data(), static_cast<std::streamsize>(filebuf.size()));
    const int n = PARAM.globalv.nlocal;
    // convert all content to string first
    // then write to ofs at once
    // this way can reduce the ofs write time significantly
    // but will consume more memory

    std::vector<std::string> row_vals(n);
    std::vector<std::string> row_cols(n);
    std::vector<int> row_nnz(n, 0);

    #pragma omp parallel for schedule(static)
    for (int row = 0; row < n; ++row)
    {
        if (pv.global2local_row(row) < 0) continue;

        auto it = XR.find(static_cast<size_t>(row));
        if (it == XR.end()) continue;

        const auto& inner = it->second;
        const int nnz = static_cast<int>(inner.size());
        row_nnz[row] = nnz;

        auto& vs = row_vals[row];
        auto& cs = row_cols[row];

        // 估算 reserve：每个 double 大约 " %.8e" 约 1+1+8+1+3+? ≈ 16~20 字符
        vs.reserve(static_cast<size_t>(nnz) * 20);
        cs.reserve(static_cast<size_t>(nnz) * 16);

        for (const auto& kv : inner)
        {
            // value（严格匹配你 write_data 的空格 + scientific(8)）
            append_sci8(vs, kv.second);

            // col index（匹配你后续 ofs << " " << i）
            cs.push_back(' ');
            char ibuf[32];
            int n = std::snprintf(ibuf, sizeof(ibuf), "%zu", kv.first);
            cs.append(ibuf, ibuf + n);
        }
    }

    // indptr
    std::vector<int> indptr(n + 1, 0);
    for (int i = 0; i < n; ++i) indptr[i + 1] = indptr[i] + row_nnz[i];

    for (int row = 0; row < n; ++row)
    {
        const auto& s = row_vals[row];
        if (!s.empty()) ofs.write(s.data(), static_cast<std::streamsize>(s.size()));
    }

    ofs.put('\n');

    for (int row = 0; row < n; ++row)
    {
        const auto& s = row_cols[row];
        if (!s.empty()) ofs.write(s.data(), static_cast<std::streamsize>(s.size()));
    }
    ofs.put('\n');

    std::string ip;
    ip.reserve(static_cast<size_t>(n + 1) * 12);
    for (int i = 0; i <= n; ++i)
    {
        ip.push_back(' ');
        char ibuf[32];
        int n = std::snprintf(ibuf, sizeof(ibuf), "%d", indptr[i]);
        ip.append(ibuf, ibuf + n);
    }
    ip.push_back('\n');
    ofs.write(ip.data(), static_cast<std::streamsize>(ip.size()));
    ModuleBase::timer::tick("ModuleIO", "output_single_R");
}

template <typename T>
void ModuleIO::output_single_R_non_reduce_binary(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, T>>& XR,
    const double& sparse_threshold,
    const Parallel_Orbitals& pv)
{
    ModuleBase::timer::tick("ModuleIO", "output_single_R");
    const int n = PARAM.globalv.nlocal;

    // Pre-compute per-row data
    std::vector<std::vector<T>> row_values(n);
    std::vector<std::vector<int>> row_cols(n);

    #pragma omp parallel for schedule(static)
    for (int row = 0; row < n; ++row)
    {
        if (pv.global2local_row(row) < 0) continue;

        auto it = XR.find(static_cast<size_t>(row));
        if (it == XR.end()) continue;

        const auto& inner = it->second;
        // Reserve space to avoid reallocations
        row_values[row].reserve(inner.size());
        row_cols[row].reserve(inner.size());

        for (const auto& kv : inner)
        {
            if (std::abs(kv.second) > sparse_threshold)
            {
                row_values[row].push_back(kv.second);
                row_cols[row].push_back(static_cast<int>(kv.first));
            }
        }
    }

    // Build indptr and flatten data
    std::vector<int> indptr(n + 1, 0);
    size_t total_nnz = 0;
    for (int row = 0; row < n; ++row)
    {
        total_nnz += row_values[row].size();
        indptr[row + 1] = indptr[row] + static_cast<int>(row_values[row].size());
    }

    // Flatten values and col_indices
    std::vector<T> values;
    std::vector<int> col_indices;
    values.reserve(total_nnz);
    col_indices.reserve(total_nnz);

    for (int row = 0; row < n; ++row)
    {
        values.insert(values.end(), row_values[row].begin(), row_values[row].end());
        col_indices.insert(col_indices.end(), row_cols[row].begin(), row_cols[row].end());
    }

    // Write in CSR format: values, col_indices, indptr
    ofs.write(reinterpret_cast<char*>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(T)));
    ofs.write(reinterpret_cast<char*>(col_indices.data()),
              static_cast<std::streamsize>(col_indices.size() * sizeof(int)));
    ofs.write(reinterpret_cast<char*>(indptr.data()),
              static_cast<std::streamsize>(indptr.size() * sizeof(int)));
    ModuleBase::timer::tick("ModuleIO", "output_single_R");
}

template void ModuleIO::output_single_R<double>(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, double>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce);

template void ModuleIO::output_single_R<std::complex<double>>(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, std::complex<double>>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce);

