#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <string>
#include <cctype>

using namespace std;

double extractCGfromSequence(const string& sequence) {
    if (sequence.empty()) return 0.0;
    int total_gc = 0;
    for (char c : sequence) {
        c = toupper(c);
        if (c == 'C' || c == 'G') total_gc++;
    }
    return static_cast<double>(total_gc) / sequence.length();
}

double extractHydrophobicAAFromORF(const string& ORFsequence) {
    if (ORFsequence.empty()) return 0.0;
    static const unordered_set<char> hydrophobic_aa = {'G', 'A', 'L', 'I', 'F', 'W', 'M', 'V', 'P'};
    int total = count_if(ORFsequence.begin(), ORFsequence.end(),
                         [](char aa) { return hydrophobic_aa.count(aa) > 0; });
    return static_cast<double>(total) / ORFsequence.length();
}

unordered_map<string, double> computeKmerFrequencies(const string& sequence, int k) {
    unordered_map<string, int> counts;
    int total_kmers = 0;

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        string kmer = sequence.substr(i, k);
        transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
        
        if (any_of(kmer.begin(), kmer.end(), [](char c){ return c != 'A' && c != 'C' && c != 'G' && c != 'T'; }))
            continue;
        
            counts[kmer]++;
        total_kmers++;
    }

    unordered_map<string, double> freqs;
    for (const auto& p : counts)
        freqs[p.first] = static_cast<double>(p.second) / total_kmers;

    return freqs;
}

unordered_map<string, vector<size_t>> findAllMotifPositions(const string& sequence, const vector<string>& motifs) {
    unordered_map<string, vector<size_t>> motif_positions;
    unordered_map<size_t, vector<string>> motifs_by_length;

    for (const string& motif : motifs)
        motifs_by_length[motif.length()].push_back(motif);

    for (size_t i = 0; i < sequence.length(); ++i) {
        for (const auto& [len, motif_list] : motifs_by_length) {
            if (i + len > sequence.length()) continue;
            string sub = sequence.substr(i, len);
            for (const string& motif : motif_list)
                if (sub == motif)
                    motif_positions[motif].push_back(i);
        }
    }

    return motif_positions;
}

pair<double, double> computeAverageAndMedianDistances(const vector<size_t>& positions) {
    if (positions.size() < 2) return {0.0, 0.0};

    vector<int> distances;
    for (size_t i = 1; i < positions.size(); ++i)
        distances.push_back(static_cast<int>(positions[i] - positions[i - 1]));

    double average = accumulate(distances.begin(), distances.end(), 0.0) / distances.size();

    sort(distances.begin(), distances.end());
    double median;

    size_t mid = distances.size() / 2;
    if (distances.size() % 2 == 0) {
        median = (distances[mid - 1] + distances[mid]) / 2.0;
    } else {
        median = distances[mid];
    }

    return {average, median};
}

vector<string> splitByComma(const string& input) {
    vector<string> result;
    stringstream ss(input); string item;
    while (getline(ss, item, ',')) {
        if (!item.empty())
            result.emplace_back(move(item));
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 5 || string(argv[1]) != "-i" || string(argv[3]) != "-t") {
        cerr << "Usage: " << argv[0] << " -i <tsv_file> -t <motif1,motif2,...>\n";
        return 1;
    }

    string tsv_file = argv[2];
    vector<string> motifs = splitByComma(argv[4]);
    vector<string> motif_features;
    for (const string& m : motifs) {
        motif_features.push_back("avgDist_" + m);
        motif_features.push_back("medianDist_" + m);
    }

    vector<string> freq_features;
    for (const string& base1 : {"A", "C", "G", "T"})
        for (const string& base2 : {"A", "C", "G", "T"})
            freq_features.push_back("freq_" + base1 + base2);

    for (const string& base1 : {"A", "C", "G", "T"})
        for (const string& base2 : {"A", "C", "G", "T"})
            for (const string& base3 : {"A", "C", "G", "T"})
                freq_features.push_back("freq_" + base1 + base2 + base3);

    // Print header
    cout << "orf_length,compression,number_hidrophobic_aa,gc_content,strand";
    for (const auto& f : motif_features) cout << "," << f;
    for (const auto& f : freq_features) cout << "," << f;
    cout << endl;

    ifstream infile(tsv_file);
    string line;
    while (getline(infile, line)) {
        stringstream ss(line);
        string sequence, longest_orf, compression, strand_s;
        getline(ss, sequence, '\t');
        getline(ss, longest_orf, '\t');
        getline(ss, compression, '\t');
        getline(ss, strand_s, '\t');

        double gc = extractCGfromSequence(sequence);
        double hydrophobic = extractHydrophobicAAFromORF(longest_orf);
        double orf_ratio = longest_orf.empty() ? 0.0 : (double)longest_orf.length()/sequence.length();
        double strand = stod(strand_s);

        auto motif_pos = findAllMotifPositions(sequence, motifs);
        unordered_map<string, pair<double, double>> motif_stats;
        for (const string& m : motifs){
            auto [avg, med] = computeAverageAndMedianDistances(motif_pos[m]);
            avg /= sequence.length();
            med /= sequence.length();
            motif_stats[m] = {avg, med};
        }
            

        auto freq2 = computeKmerFrequencies(sequence, 2);
        auto freq3 = computeKmerFrequencies(sequence, 3);

        cout << orf_ratio << "," << compression << "," << hydrophobic << "," << gc << "," << strand;
        for (const string& m : motifs)
            cout << "," << motif_stats[m].first << "," << motif_stats[m].second;

        for (const string& base1 : {"A", "C", "G", "T"}) {
            for (const string& base2 : {"A", "C", "G", "T"}) {
                string kmer = base1 + base2;
                cout << "," << (freq2.count(kmer) ? freq2.at(kmer) : 0.0);
            }
        }

        for (const string& base1 : {"A", "C", "G", "T"}) {
            for (const string& base2 : {"A", "C", "G", "T"}) {
                for (const string& base3 : {"A", "C", "G", "T"}) {
                    string kmer = base1 + base2 + base3;
                    cout << "," << (freq3.count(kmer) ? freq3.at(kmer) : 0.0);
                }
            }
        }

        cout << endl;
    }
    return 0;
}
