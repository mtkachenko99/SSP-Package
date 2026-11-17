#pragma once
// Запись результатов в CSV (заголовок-комментарий с единицами и происхождением).

#include <string>
#include <vector>

namespace ssp::app {

class CsvWriter {
 public:
  CsvWriter(const std::string& path, const std::vector<std::string>& comments, const std::vector<std::string>& columns);
  ~CsvWriter();
  void row(const std::vector<double>& values);
  void row_text(const std::vector<std::string>& values);

 private:
  void* f_ = nullptr;
};

void make_dirs(const std::string& path);
std::string join_path(const std::string& a, const std::string& b);

}  // namespace ssp::app
