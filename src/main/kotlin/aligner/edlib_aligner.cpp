#include <cstdio>
#include "edlib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {

    for (int i = 1; i < argc; i++) {
        string dir = argv[i];
        cout << "Aligning " << dir << "...";
        cout.flush();

        ifstream ref(dir + "/_reference");
        string reference, part;
        while (ref >> part) {
            reference += part;
        }
        ref.close();

        ifstream reads(dir + "/_reads");
        ofstream align(dir + "/_aligned");
        while (reads >> part) {
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

            align << result << ' ' << start << '\n';
        }
        reads.close();
        align.close();

        cout << " done" << endl;
        cout.flush();
    }

}