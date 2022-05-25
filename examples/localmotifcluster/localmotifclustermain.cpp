// localmotifclustermain.cpp : Defines the entry point for the console application.
//
#include "localmotifcluster.h"
#include "stdafx.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

// Hash function
struct hashFunction {
    size_t operator()(const vector<int>& myVector) const
    {
        std::hash<int> hasher;
        size_t answer = 0;

        for (int i : myVector) {
            answer ^= hasher(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
        }
        return answer;
    }
};

template <class BidiIter>
BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random)
{
    size_t left = std::distance(begin, end);
    while (num_random--) {
        BidiIter r = begin;
        std::advance(r, rand() % left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

vector<unordered_set<int>> read_community(string path)
{
    ifstream file(path);
    if (!file.good()) {
        cerr << "file is not good" << endl;
        exit(0);
    }
    vector<unordered_set<int>> communities;
    unordered_set<vector<int>, hashFunction> nodups;
    string line;
    while (getline(file, line)) {
        unordered_set<int> community;
        istringstream ss(line);
        int num;
        vector<int> cmvec;
        while (ss >> num) {
            community.insert(num);
            cmvec.push_back(num);
        }
        sort(cmvec.begin(), cmvec.end());
        if (nodups.find(cmvec) == nodups.end()) {
            communities.push_back(community);
            nodups.insert(cmvec);
        }
    }
    return communities;
}

string convert_from_tstr(const TStr tstr)
{
    string s;
    int len = tstr.Len();
    for (int i = 0; i < tstr.Len(); ++i) {
        s += tstr.GetCh(i);
    }
    return s;
}

int main(int argc, char* argv[])
{
    auto start_total = high_resolution_clock::now();
    auto start_motif = high_resolution_clock::now();
    Env = TEnv(argc, argv, TNotify::StdNotify);
    Env.PrepArgs(TStr::Fmt("Local motif clustering. build: %s, %s. Time: %s",
        __TIME__, __DATE__, TExeTm::GetCurTm()));
    TExeTm ExeTm;

    Try

        const TStr graph_filename
        = Env.GetIfArgPrefixStr("-i:", "graph/karate_club.gr", "Input graph file");
    const bool IsDirected = Env.GetIfArgPrefixBool("-d:", false, "Directed graph?");

    // community file
    const TStr community_filename = Env.GetIfArgPrefixStr("-c:", "community/karate.cmty.txt", "Input community file");
    string community_path = convert_from_tstr(community_filename);
    vector<unordered_set<int>> communities = read_community(community_path);

    // output file
    const TStr output_filename = Env.GetIfArgPrefixStr("-o:", "output/karate.output.txt", "Output file");
    string output_path = convert_from_tstr(output_filename);
    ofstream output(output_path);

    ProcessedGraph graph_p;
    if (IsDirected) {
        const TStr motif = Env.GetIfArgPrefixStr("-m:", "triad", "Motif type");
        MotifType mt = ParseMotifType(motif, IsDirected);
        PNGraph graph;
        if (graph_filename.GetFExt().GetLc() == ".ngraph") {
            TFIn FIn(graph_filename);
            graph = TNGraph::Load(FIn);
        } else if (graph_filename.GetFExt().GetLc() == ".ungraph") {
            TExcept::Throw("Warning: input graph is an undirected graph!!");
        } else {
            graph = TSnap::LoadEdgeList<PNGraph>(graph_filename, 0, 1);
        }
        TSnap::DelSelfEdges(graph);
        graph_p = ProcessedGraph(graph, mt);
    } else {
        const TStr motif = Env.GetIfArgPrefixStr("-m:", "clique3", "Motif type");
        MotifType mt = ParseMotifType(motif, IsDirected);
        PUNGraph graph;
        if (graph_filename.GetFExt().GetLc() == ".ungraph") {
            TFIn FIn(graph_filename);
            graph = TUNGraph::Load(FIn);
        } else if (graph_filename.GetFExt().GetLc() == ".ngraph") {
            TExcept::Throw("Warning: input graph is a directed graph!!");
        } else {
            graph = TSnap::LoadEdgeList<PUNGraph>(graph_filename, 0, 1);
        }
        TSnap::DelSelfEdges(graph);
        graph_p = ProcessedGraph(graph, mt);
    }
    const TFlt alpha = Env.GetIfArgPrefixFlt("-a:", 0.98, "alpha");
    const TFlt eps = Env.GetIfArgPrefixFlt("-e:", 0.0001, "eps");
    // const TInt seed = Env.GetIfArgPrefixInt("-s:", 1, "Seed");
    // mappr.computeAPPR(graph_p, seed, alpha, eps / graph_p.getTotalVolume() * graph_p.getTransformedGraph()->GetNodes());
    // mappr.sweepAPPR(-1);
    // mappr.printProfile();
    auto end_motif = high_resolution_clock::now();
    std::cout << "motif discovery time: " << (double)duration_cast<microseconds>(end_motif - start_motif).count() / 1000000 << endl;

    output << "id,precision,recall,f1" << endl;
    double sum_precision = 0, sum_recall = 0, sum_f1 = 0;
    unordered_map<int, vector<int>> nd2cluster;
    for (int i = 0; i < communities.size(); ++i) {
        auto community = communities[i];
        double max_p = 0, max_r = 0, max_f1 = 0;
        int best_seed = -1;
        for (auto seed : community) {
            vector<int> expectations;
            if (nd2cluster.find(seed) != nd2cluster.end())
                expectations = nd2cluster[seed];
            else {
                MAPPR mappr;
                mappr.computeAPPR(graph_p, seed, alpha, eps / graph_p.getTotalVolume() * graph_p.getTransformedGraph()->GetNodes());
                mappr.sweepAPPR(-1);
                expectations = mappr.getNodesInOrder(graph_p, true);
                nd2cluster[seed] = expectations;
            }
            int hit = 0;
            for (auto nd : expectations) {
                if (community.find(nd) != community.end())
                    ++hit;
            }
            double precision = (double)hit / expectations.size();
            double recall = (double)hit / community.size();
            double f1 = 2 * precision * recall / (precision + recall);
            if (f1 > max_f1) {
                max_p = precision;
                max_r = recall;
                max_f1 = f1;
                best_seed = seed;
            }
        }
        output << i << "," << max_p << "," << max_r << "," << max_f1 << endl;
        std::cout << i << "," << max_p << "," << max_r << "," << max_f1 << endl;
        sum_precision += max_p;
        sum_recall += max_r;
        sum_f1 += max_f1;
    }
    int len = communities.size();
    std::cout << "FINAL precision: " << sum_precision / len << ", recall: " << sum_recall / len << ", f1: " << sum_f1 / len << endl;

    Catch
        printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(),
            TSecTm::GetCurTm().GetTmStr().CStr());

    auto end_total = high_resolution_clock::now();
    std::cout << "total time: " << (double)duration_cast<microseconds>(end_total - start_total).count() / 1000000 << endl;

    return 0;
}
