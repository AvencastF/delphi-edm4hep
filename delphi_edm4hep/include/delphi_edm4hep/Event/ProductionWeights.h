#pragma once

#include <cmath>
#include <cstdlib>
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
  // Optional generator tau polarimeters (prepare_weights.py TAU_POLARIMETERS section), per event.
  struct Taus { std::vector<int> pdg, meMode; std::vector<double> p4, polarimeter; };
  int run_{};
  double sigma_{}, error_{};
  std::string provenance_;
  std::map<int, Entry> entries_;
  std::set<int> converted_;
  bool hasTaus_{};
  std::map<int, Taus> taus_;
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
    if (in >> extra && extra == "TAU_POLARIMETERS") {
      int n{};
      if (!(in >> n) || n < 0) throw std::runtime_error("Invalid tau polarimeter header");
      hasTaus_ = true;
      for (int i = 0; i < n; ++i) {
        int event{}, pdg{}, meMode{}; double v[8];
        if (!(in >> event >> pdg >> meMode) || !entries_.count(event) || std::abs(pdg) != 15)
          throw std::runtime_error("Invalid tau polarimeter record");
        // Tokens via stod: an unidentified channel (meMode -1) carries nan, which operator>> rejects.
        for (double& x : v) {
          std::string token;
          if (!(in >> token)) throw std::runtime_error("Invalid tau polarimeter record");
          x = std::stod(token);
        }
        Taus& t = taus_[event];
        t.pdg.push_back(pdg); t.meMode.push_back(meMode);
        t.p4.insert(t.p4.end(), v, v + 4); t.polarimeter.insert(t.polarimeter.end(), v + 4, v + 8);
      }
      if (in >> extra) throw std::runtime_error("Extra tau polarimeter records");
    } else if (in) throw std::runtime_error("Extra production weight records");
  }

  template<class Frame> void event(Frame& frame, int run, int event) {
    const auto it = entries_.find(event);
    if (run != run_ || it == entries_.end() || !converted_.insert(event).second)
      throw std::runtime_error("Missing, wrong-run or duplicate MC weight identity");
    frame.putParameter("mc_gen_weight", it->second.weight);
    frame.putParameter("mc_generator_trial", it->second.trial);
    if (!hasTaus_) return;
    // Every event of a chunk with the section carries the four parallel lists, empty without decayed taus.
    const Taus& t = taus_[event];
    frame.putParameter("mc_tau_pdg", t.pdg);
    frame.putParameter("mc_tau_me_mode", t.meMode);
    frame.putParameter("mc_tau_p4", t.p4);
    frame.putParameter("mc_tau_polarimeter", t.polarimeter);
  }

  template<class Frame> void metadata(Frame& frame) const {
    frame.putParameter("mc_weight_schema", 1);
    frame.putParameter("mc_weight_names", std::vector<std::string>{"nominal"});
    frame.putParameter("mc_run_number", run_);
    frame.putParameter("mc_cross_section_pb", sigma_);
    frame.putParameter("mc_cross_section_error_pb", error_);
    frame.putParameter("mc_production_json", provenance_);
    if (hasTaus_) frame.putParameter("mc_tau_polarimeter_schema", 1);
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
