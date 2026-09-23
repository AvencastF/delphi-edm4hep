#pragma once

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace delphi_edm4hep {
// No global state or dependency on the legacy commons. Frame is templated so
// identity and signed-weight bookkeeping can also be tested without ROOT.
class ProductionWeights {
  struct Entry { int trial; double weight; };
  int run_{};
  double sigma_{}, error_{};
  std::string provenance_;
  std::map<int, Entry> entries_;
  std::set<int> converted_;
public:
  explicit ProductionWeights(const std::string& path) {
    std::ifstream in(path);
    std::string magic;
    int count{};
    if (!std::getline(in, magic) || magic != "DELPHI_MC_WEIGHTS_V1" ||
        !(in >> run_ >> count >> sigma_ >> error_) || run_ >= 0 || count <= 0 ||
        !std::isfinite(sigma_) || sigma_ <= 0 || !std::isfinite(error_) || error_ < 0)
      throw std::runtime_error("Invalid production weight header");
    std::getline(in, provenance_); // numeric-header newline
    if (!std::getline(in, provenance_) || provenance_.empty())
      throw std::runtime_error("Missing production weight provenance");
    std::set<int> trials;
    for (int i = 1; i <= count; ++i) {
      int event{}, trial{}; double weight{};
      if (!(in >> event >> trial >> weight) || event != i || !std::isfinite(weight) ||
          (trial != -1 && (trial <= 0 || !trials.insert(trial).second)))
        throw std::runtime_error("Invalid production weight record");
      entries_.emplace(event, Entry{trial, weight});
    }
    std::string extra;
    if (in >> extra) throw std::runtime_error("Extra production weight records");
  }

  template<class Frame> void event(Frame& frame, int run, int event) {
    const auto it = entries_.find(event);
    if (run != run_ || it == entries_.end() || !converted_.insert(event).second)
      throw std::runtime_error("Missing, wrong-run or duplicate MC weight identity");
    frame.putParameter("mc_gen_weight", it->second.weight);
    frame.putParameter("mc_generator_trial", it->second.trial);
  }

  template<class Frame> void metadata(Frame& frame) const {
    frame.putParameter("mc_weight_schema", 1);
    frame.putParameter("mc_weight_names", std::vector<std::string>{"nominal"});
    frame.putParameter("mc_run_number", run_);
    frame.putParameter("mc_cross_section_pb", sigma_);
    frame.putParameter("mc_cross_section_error_pb", error_);
    frame.putParameter("mc_production_json", provenance_);
    std::vector<int> generated, trials, converted, missing, reserve;
    std::vector<double> weights;
    const int last = converted_.empty() ? 0 : *converted_.rbegin();
    for (const auto& [id, entry] : entries_) {
      generated.push_back(id); trials.push_back(entry.trial); weights.push_back(entry.weight);
      if (converted_.count(id)) converted.push_back(id);
      else if (id <= last) missing.push_back(id);
      else reserve.push_back(id);
    }
    frame.putParameter("mc_generated_event_ids", generated);
    frame.putParameter("mc_generated_trial_ids", trials);
    frame.putParameter("mc_generated_weights", weights);
    frame.putParameter("mc_converted_event_ids", converted);
    frame.putParameter("mc_missing_ids_before_last_output", missing);
    frame.putParameter("mc_ids_after_last_output", reserve);
    // Last-output accounting is an observation, not proof that later records
    // were never attempted. Keep this distinct from a simulation efficiency.
    auto sums = [&](const std::string& prefix, const std::vector<int>& ids) {
      double sum = 0, square = 0, absolute = 0; int negative = 0;
      for (int id : ids) {
        const double w = entries_.at(id).weight;
        sum += w; square += w*w; absolute += std::abs(w); negative += w < 0;
      }
      if (!std::isfinite(sum) || !std::isfinite(square) || !std::isfinite(absolute))
        throw std::runtime_error("Nonfinite production weight sums");
      frame.putParameter(prefix + "count", static_cast<int>(ids.size()));
      frame.putParameter(prefix + "sumw", sum);
      frame.putParameter(prefix + "sumw2", square);
      frame.putParameter(prefix + "sumabsw", absolute);
      frame.putParameter(prefix + "negative_count", negative);
    };
    sums("mc_generated_", generated);
    sums("mc_converted_", converted);
    sums("mc_missing_", missing);
    sums("mc_after_last_output_", reserve);
  }
};
} // namespace delphi_edm4hep
