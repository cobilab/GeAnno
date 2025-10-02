#include <iostream>
#include <vector>
#include <cctype>
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstring>
#include <numeric>
#include <unordered_map>

using namespace std;

double extractCGfromSequence(const string& sequence) {
    if (sequence.empty()) return 0.0;
    int total_gc = count_if(sequence.begin(), sequence.end(), [](char c) {
        c = toupper(c);
        return c == 'C' || c == 'G';
    });
    return static_cast<double>(total_gc) / sequence.length();
}

unordered_map<string, double> computeKmerFrequencies(const string& sequence, int k) {
    unordered_map<string, int> counts;
    int total_kmers = 0;

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        string kmer = sequence.substr(i, k);
        transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);

        // Skip non-ACGT characters
        if (any_of(kmer.begin(), kmer.end(), [](char c){ return c != 'A' && c != 'C' && c != 'G' && c != 'T'; }))
            continue;

        counts[kmer]++;
        total_kmers++;
    }

    unordered_map<string, double> frequencies;
    for (const auto& pair : counts) {
        frequencies[pair.first] = static_cast<double>(pair.second) / total_kmers;
    }

    return frequencies;
}

double extractHydrophobicAAFromORF(const string& ORFsequence) {
    if (ORFsequence.empty()) return 0.0;
    static const unordered_set<char> hydrophobic_aa = {'G', 'A', 'L', 'I', 'F', 'W', 'M', 'V', 'P'};
    int total = count_if(ORFsequence.begin(), ORFsequence.end(),
                         [](char aa) { return hydrophobic_aa.count(aa) > 0; });
    return static_cast<double>(total) / ORFsequence.length();
}

unordered_map<string, vector<size_t>> findAllMotifPositions(const string& sequence, const vector<string>& motifs) {
    unordered_map<string, vector<size_t>> motif_positions;
    for (const string& motif : motifs) {
        for (size_t pos = 0; pos + motif.size() <= sequence.size(); ++pos) {
            if (sequence.substr(pos, motif.size()) == motif)
                motif_positions[motif].push_back(pos);
        }
    }
    return motif_positions;
}

pair<double, double> computeAverageAndMedianDistances(const vector<size_t>& positions) {
    if (positions.size() < 2) return {0.0, 0.0};
    vector<int> distances;
    for (size_t i = 1; i < positions.size(); ++i)
        distances.push_back(positions[i] - positions[i - 1]);

    double avg = accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
    sort(distances.begin(), distances.end());
    double median = distances.size() % 2 == 0
        ? (distances[distances.size()/2 - 1] + distances[distances.size()/2]) / 2.0
        : distances[distances.size()/2];
    return {avg, median};
}

void initializeFile(ofstream& out, const vector<string>& motifs) {
    out << "species,label,orf_length,compression,number_hidrophobic_aa,gc_content";
    
    // dist motifs
    for (const string& motif : motifs)
        out << ",avgDist_" << motif << ",medianDist_" << motif;

    // 2-kmer e 3-kmer frequencies
    for (const string& base1 : {"A", "C", "G", "T"})
        for (const string& base2 : {"A", "C", "G", "T"})
            out << ",freq_" << base1 << base2;

    for (const string& base1 : {"A", "C", "G", "T"})
        for (const string& base2 : {"A", "C", "G", "T"})
            for (const string& base3 : {"A", "C", "G", "T"})
                out << ",freq_" << base1 << base2 << base3;

    out << "\n"; 

}

void writeToFile(ofstream& out,
    const string& species,
    const string& label,
    const double& orf_ratio,
    const string& compression,
    const double& hydro_ratio,
    const double& gc_ratio,
    const vector<pair<double, double>>& motif_stats,
    const unordered_map<string, double>& kmer2,
    const unordered_map<string, double>& kmer3
    ) {

    out << species << "," << label << "," << orf_ratio << "," << compression << ","<< hydro_ratio << "," << gc_ratio;
    
    for (const auto& [avg, med] : motif_stats)
        out << "," << avg << "," << med;

    for (const string& base1 : {"A", "C", "G", "T"}) {
        for (const string& base2 : {"A", "C", "G", "T"}) {
            string kmer = base1 + base2;
            out << "," << (kmer2.count(kmer) ? kmer2.at(kmer) : 0.0);
        }
    }

    for (const string& base1 : {"A", "C", "G", "T"}) {
        for (const string& base2 : {"A", "C", "G", "T"}) {
            for (const string& base3 : {"A", "C", "G", "T"}) {
                string kmer = base1 + base2 + base3;
                out << "," << (kmer3.count(kmer) ? kmer3.at(kmer) : 0.0);
            }
        }
    }
    
    out << "\n";
}

void iterateThroughCSV(ofstream& out, const string& filename, const vector<string>& motifs) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open input CSV.\n";
        return;
    }

    string line;
    getline(file, line);  // Skip header

    while (getline(file, line)) {
        vector<string> fields;
        stringstream ss(line);
        string item;
        while (getline(ss, item, ',')) fields.push_back(item);

        //"species","label","seq","longest_seq","compression_nc"
        string species = fields[0];
        string label = fields[1];
        string seq = fields[2];
        string orf = fields[3];
        string compression = fields[4];
        //string strand = fields[5];
    
        double gc = extractCGfromSequence(seq);
        double hydro = extractHydrophobicAAFromORF(orf);
        double orf_ratio = static_cast<double>(orf.length()) / seq.length();

        auto motif_positions = findAllMotifPositions(seq, motifs);
        vector<pair<double, double>> motif_stats;
        for (const string& motif : motifs) {
            auto [avg, med] = computeAverageAndMedianDistances(motif_positions[motif]);
            motif_stats.emplace_back(avg / seq.length(), med / seq.length());
        }

        auto kmer2 = computeKmerFrequencies(seq, 2);
        auto kmer3 = computeKmerFrequencies(seq, 3);

        writeToFile(out, species, label, orf_ratio, compression, hydro, gc, motif_stats, kmer2, kmer3);
    }

    file.close();
}

vector<string> splitByComma(const string& input) {
    vector<string> result;
    stringstream ss(input);
    string token;
    while (getline(ss, token, ',')) {
        if (!token.empty()) result.push_back(token);
    }
    return result;
}

int main(int argc, char* argv[]) {
    string input_csv = "", output_csv = "";
    vector<string> motifs;

    for (int i = 1; i < argc - 1; ++i) {
        if (strcmp(argv[i], "-i") == 0) input_csv = argv[i + 1];
        else if (strcmp(argv[i], "-o") == 0) output_csv = argv[i + 1];
        else if (strcmp(argv[i], "-t") == 0) motifs = splitByComma(argv[i + 1]);
    }

    ofstream out(output_csv);
    if (!out.is_open()) {
        cerr << "Could not open output file.\n";
        return 1;
    }

    initializeFile(out, motifs);
    iterateThroughCSV(out, input_csv, motifs);
    out.close();
    return 0;
}
