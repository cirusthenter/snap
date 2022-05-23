// localmotifclustermain.cpp : Defines the entry point for the console application.
//
#include "localmotifcluster.h"
#include "stdafx.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

vector<unordered_set<int>> read_community(string path)
{
    ifstream file(path);
    if (!file.good()) {
        cerr << "file is not good" << endl;
        exit(0);
    }
    vector<unordered_set<int>> communities;
    string line;
    while (getline(file, line)) {
        unordered_set<int> community;
        istringstream ss(line);
        int num;
        while (ss >> num)
            community.insert(num);
        communities.push_back(community);
    }
    return communities;
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

        const bool IsDirected
        = Env.GetIfArgPrefixBool("-d:", false, "Directed graph?");

    ProcessedGraph graph_p;
    if (IsDirected) {
        const TStr graph_filename = Env.GetIfArgPrefixStr("-i:", "C-elegans-frontal.txt", "Input graph file");
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
        const TStr graph_filename = Env.GetIfArgPrefixStr("-i:", "C-elegans-frontal.txt", "Input graph file");
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
        cout << "#nodes: " << graph->GetNodes() << ", "
             << "#edges: " << graph->GetEdges() << endl;
        graph_p = ProcessedGraph(graph, mt);
    }
    const TFlt alpha = Env.GetIfArgPrefixFlt("-a:", 0.98, "alpha");
    const TFlt eps = Env.GetIfArgPrefixFlt("-e:", 0.0001, "eps");
    // const TInt seed = Env.GetIfArgPrefixInt("-s:", 1, "Seed");
    // mappr.computeAPPR(graph_p, seed, alpha, eps / graph_p.getTotalVolume() * graph_p.getTransformedGraph()->GetNodes());
    // mappr.sweepAPPR(-1);
    // mappr.printProfile();
    auto end_motif = high_resolution_clock::now();
    cout << "motif discovery time: " << (double)duration_cast<microseconds>(end_motif - start_motif).count() / 1000000 << endl;

    vector<unordered_set<int>> communities = read_community("community/karate.cmty.txt");
    for (auto community : communities) {
        double max_p = 0, max_r = 0, max_f1 = 0;
        for (auto seed : community) {
            MAPPR mappr;
            mappr.computeAPPR(graph_p, seed, alpha, eps / graph_p.getTotalVolume() * graph_p.getTransformedGraph()->GetNodes());
            mappr.sweepAPPR(-1);
            vector<int> expectations = mappr.getNodesInOrder();
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
            }
        }
        cout << "precision: " << max_p << ", recall: " << max_r << ", f1: " << max_f1 << endl;
    }

    Catch
        printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(),
            TSecTm::GetCurTm().GetTmStr().CStr());

    auto end_total = high_resolution_clock::now();
    cout << "total time: " << (double)duration_cast<microseconds>(end_total - start_total).count() / 1000000 << endl;

    return 0;
}
