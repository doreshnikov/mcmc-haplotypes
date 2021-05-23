#include <cstdio>
#include "edlib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include "util/ThreadPool.h"

using namespace std;

int main(int argc, char *argv[]) {

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;

    int nthreads = thread::hardware_concurrency();
    ThreadPool tpool(nthreads);

    for (int i = 1; i < argc; i++) {
        string dir = argv[i];
        cout << "Aligning " << dir << "...";
        cout.flush();
        auto t1 = high_resolution_clock::now();

        ifstream ref(dir + "/_reference");
        string reference, part;
        while (ref >> part) {
            reference += part;
        }
        ref.close();

        ifstream reads(dir + "/_reads");
        vector<future<pair<string, int>>> promises;
        while (reads >> part) {
            promises.emplace_back(
                tpool.enqueue([part, &reference]() {
                    string rev = part;
                    reverse(rev.begin(), rev.end());
                    string cmp = part;
                    for (char &idx : cmp) {
                        if (idx == 'A') idx = 'T';
                        else if (idx == 'T') idx = 'A';
                        else if (idx == 'C') idx = 'G';
                        else if (idx == 'G') idx = 'C';
                    }
                    string revcmp = cmp;
                    reverse(revcmp.begin(), revcmp.end());

                    vector<string> opts{part, rev, cmp, revcmp};
                    string result;
                    int distance = static_cast<int>(part.length()), start = 0;
                    for (string const &opt : opts) {
                        EdlibAlignResult res = edlibAlign(
                                opt.c_str(), static_cast<int>(opt.length()),
                                reference.c_str(), static_cast<int>(reference.length()),
                                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0)
                        );
                        if (res.editDistance < distance) {
                            distance = res.editDistance;
                            start = res.startLocations[0];
                            result = opt;
                        }
                        edlibFreeAlignResult(res);
                    }
                    return make_pair(result, start);
                })
            );
        }
        reads.close();

        ofstream align(dir + "/_aligned");
        for (auto &promise : promises) {
            auto result = promise.get();
            align << result.first << ' ' << result.second << '\n';
        }
        align.close();

        auto t2 = high_resolution_clock::now();
        auto time_ms = duration_cast<milliseconds>(t2 - t1);
        cout << " done in " << time_ms.count() / 1000 << "s\n";
        cout.flush();
    }

}