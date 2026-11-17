#include "output.hpp"

#include <cstdio>
#include <filesystem>
#include <stdexcept>

namespace ssp::app {

CsvWriter::CsvWriter(const std::string& path, const std::vector<std::string>& comments,
                     const std::vector<std::string>& columns) {
  FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw std::runtime_error("cannot create file: " + path);
  f_ = f;
  for (const auto& c : comments) std::fprintf(f, "# %s\n", c.c_str());
  for (std::size_t i = 0; i < columns.size(); ++i) std::fprintf(f, i ? ",%s" : "%s", columns[i].c_str());
  std::fprintf(f, "\n");
}

CsvWriter::~CsvWriter() {
  if (f_) std::fclose(static_cast<FILE*>(f_));
}

void CsvWriter::row(const std::vector<double>& v) {
  FILE* f = static_cast<FILE*>(f_);
  for (std::size_t i = 0; i < v.size(); ++i) std::fprintf(f, i ? ",%.10g" : "%.10g", v[i]);
  std::fprintf(f, "\n");
}

void CsvWriter::row_text(const std::vector<std::string>& v) {
  FILE* f = static_cast<FILE*>(f_);
  for (std::size_t i = 0; i < v.size(); ++i) std::fprintf(f, i ? ",%s" : "%s", v[i].c_str());
  std::fprintf(f, "\n");
}

void make_dirs(const std::string& path) { std::filesystem::create_directories(path); }

std::string join_path(const std::string& a, const std::string& b) { return (std::filesystem::path(a) / b).string(); }

}  // namespace ssp::app
