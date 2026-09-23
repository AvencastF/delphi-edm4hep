import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('prepare_weights', ROOT/'scripts/production/prepare_weights.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class WeightsTest(unittest.TestCase):
    def fixture(self, path):
        (path/'task.json').write_text(json.dumps(dict(run=123, generated_requested=4, simulated_requested=2, generator='bhwide', seed=42)))
        (path/'generation.json').write_text(json.dumps(dict(generated=4, sumw=5, sumw2=15, sigma_mb=2e-6, sigma_error_mb=1e-8)))
        (path/'event-weights.csv').write_text('hepmc_event,weight,generator_trial\n1,1,3\n2,-1,5\n3,3,7\n4,2,9\n')

    def test_summary_and_identity_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            self.fixture(path)
            module.prepare(path)
            text = (path/'event-weights.csv').read_text()
            for bad in [text.replace('2,-1,5', '1,-1,5'), text.replace('2,-1,5', '2,nan,5'), text.replace('2,-1,5', '2,1,5'), text.replace('2,-1,5', '2,-1,3')]:
                (path/'event-weights.csv').write_text(bad)
                with self.assertRaises(ValueError): module.prepare(path)

    def test_cpp_signed_weights_gaps_and_frame_types(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            self.fixture(path)
            module.prepare(path)
            source = r'''
#include "delphi_edm4hep/Event/ProductionWeights.h"
#include <any>
#include <cassert>
struct Frame {
 std::map<std::string,std::any> values;
 template<class T> void putParameter(const std::string& key, T value) { values[key]=value; }
 template<class T> T get(const std::string& key) { return std::any_cast<T>(values.at(key)); }
};
int main(int argc, char** argv) {
 (void)argc;
 delphi_edm4hep::ProductionWeights weights(argv[1]);
 Frame first, second, meta;
 weights.event(first,-123,2); weights.event(second,-123,3);
 assert(first.get<double>("mc_gen_weight")==-1);
 assert(first.get<int>("mc_generator_trial")==5);
 weights.metadata(meta);
 assert(meta.get<double>("mc_generated_sumw")==5);
 assert(meta.get<double>("mc_generated_sumw2")==15);
 assert(meta.get<double>("mc_converted_sumw")==2);
 assert(meta.get<double>("mc_converted_sumw2")==10);
 assert(meta.get<int>("mc_converted_negative_count")==1);
 assert(meta.get<std::vector<int>>("mc_missing_ids_before_last_output")==std::vector<int>{1});
 assert(meta.get<std::vector<int>>("mc_ids_after_last_output")==std::vector<int>{4});
 for (auto pair : std::vector<std::pair<int,int>>{{-124,1},{-123,5},{-123,2}}) {
  bool failed=false;
  try { weights.event(first,pair.first,pair.second); } catch(const std::runtime_error&) { failed=true; }
  assert(failed);
 }
}
'''
            (path/'test.cpp').write_text(source)
            subprocess.run(['c++','-std=c++17','-Wall','-Wextra','-I'+str(ROOT/'delphi_edm4hep/include'),str(path/'test.cpp'),'-o',str(path/'test')],check=True)
            subprocess.run([str(path/'test'),str(path/'mc-weights.txt')],check=True)


if __name__ == '__main__': unittest.main()
