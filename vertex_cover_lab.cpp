#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Graph {
public:
    using Vertex = size_t;
    using VertexSet = std::unordered_set<Vertex>;
    using AdjacencyList = std::unordered_map<Vertex, VertexSet>;
    
    void AddVertex(Vertex v) {
        adjacency_list_[v];
    }
    
    void AddEdge(Vertex u, Vertex v) {
        adjacency_list_[u].insert(v);
        adjacency_list_[v].insert(u);
    }
    
    const VertexSet& AdjacentVertices(Vertex v) const {
        const auto it = adjacency_list_.find(v);
        if (it != adjacency_list_.end()) {
            return it->second;
        } else {
            return empty_set_;
        }
    }
    
    VertexSet AllVertices() const {
        VertexSet vs;
        vs.reserve(adjacency_list_.size());
        for (const auto& pair : adjacency_list_) {
            const auto& vertex = pair.first;
            vs.insert(vertex);
        }
        return vs;
    }
    
    const AdjacencyList& AsAdjacencyList() const {
        return adjacency_list_;
    }
    
private:
    AdjacencyList adjacency_list_;
    static const VertexSet empty_set_;
};

const Graph::VertexSet Graph::empty_set_;

class VertexCover {
public:
    explicit VertexCover(const Graph& graph)
    : graph_(graph), set_(graph.AllVertices()) {}
    
    Graph::VertexSet CandidatesToRemove() const {
        Graph::VertexSet candidates;
        for (const auto v : set_) {
            if (IsDeletable(v)) {
                candidates.insert(v);
            }
        }
        return candidates;
    }
    
    const Graph::VertexSet& CandidatesToAdd() const {
        return set_complement_;
    }
    
    void Add(Graph::Vertex v) {
        set_.insert(v);
        set_complement_.erase(v);
    }
    
    void Remove(Graph::Vertex v) {
        set_.erase(v);
        set_complement_.insert(v);
    }
    
    size_t Cost() const {
        return set_.size();
    }
    
    Graph::VertexSet::const_iterator find(Graph::Vertex v) const {
        return set_.find(v);
    }
    
    Graph::VertexSet::const_iterator begin() const {
        return set_.begin();
    }
    
    Graph::VertexSet::const_iterator end() const {
        return set_.end();
    }
    
    const Graph& GetGraph() const {
        return graph_;
    }
    
private:
    bool IsDeletable(Graph::Vertex v) const {
        for (const auto u : graph_.AdjacentVertices(v)) {
            if (find(u) == end()) {
                return false;
            }
        }
        return true;
    }
    
    const Graph& graph_;
    Graph::VertexSet set_;
    Graph::VertexSet set_complement_;
};

void GraphEdges(std::ostream& out, const Graph::AdjacencyList& adjacency_list) {
    for (const auto& pair : adjacency_list) {
        const auto& vertex = pair.first;
        const auto& neighbours = pair.second;
        for (const auto adj_vertex : neighbours) {
            out << "\t" << vertex << " -- " << adj_vertex << "\n";
        }
    }
}

// Use http://www.webgraphviz.com to take a look at the graph
void GraphViz(std::ostream& out, const Graph& graph) {
    out << "strict graph {\n";
    for (const auto& pair : graph.AsAdjacencyList()) {
        const auto& vertex = pair.first;
        out << "\t" << vertex << "\n";
    }
    GraphEdges(out, graph.AsAdjacencyList());
    out << "}\n";
}

void GraphViz(std::ostream& out, const VertexCover& vertex_cover) {
    out << "strict graph {\n";
    for (const auto& pair : vertex_cover.GetGraph().AsAdjacencyList()) {
        const auto& vertex = pair.first;
        if (vertex_cover.find(vertex) != vertex_cover.end()) {
            out << "\t" <<  vertex << " [shape=doublecircle]\n";
        } else {
            out << "\t" << vertex << "\n";
        }
    }
    GraphEdges(out, vertex_cover.GetGraph().AsAdjacencyList());
    out << "}\n";
}

struct DebugInfo {
    std::vector<size_t> costs;
};

// Use http://gnuplot.respawned.com/ to plot costs
std::ostream& operator<<(std::ostream& out, const DebugInfo& debug_info) {
    for (size_t i = 0; i < debug_info.costs.size(); ++i) {
        out << i << " " << debug_info.costs[i] << "\n";
    }
    return out;
}

class VertexCoverSolver {
public:
    virtual VertexCover Solve(const Graph& graph,
                              DebugInfo& debug_info) const = 0;
    virtual ~VertexCoverSolver() = default;
};

int InitRandSeed(int argc, const char* argv[]) {
    int rand_seed;
    if (argc >= 2) {
        rand_seed = atoi(argv[1]);
    } else {
        rand_seed = time(nullptr);
    }
    srand(rand_seed);
    return rand_seed;
}

class GradientDescent final : public VertexCoverSolver {
    
    bool checker(const VertexCover& vertexset, Graph::Vertex& elem_to_check) const {
        bool flag = false;
        for (auto elem : vertexset) {
            if (elem == elem_to_check) {
                flag == true;
                break;
            }
        }
        
        return flag;
    }
    
    VertexCover Solve(const Graph& graph, DebugInfo& debug_info) const {
        
        int rand_vertex = rand() % graph.AllVertices().size();
        VertexCover S(graph);
        
        auto vertex = graph.AllVertices().begin();
        
        int i = 0;
        
        for (auto elem : graph.AllVertices()){
            if (i == rand_vertex) {
                S.Add(elem);
            }
            i++;
        }
        
        while (true) {
            
            Graph::VertexSet TempVertexesAddition;
            
            for (auto vert = S.begin(); vert != S.end(); ++vert) {
                Graph::VertexSet temp = graph.AdjacentVertices(*vert);
                
                for (auto elem : temp) {
                    
                    if (checker(S, elem)) {
                        TempVertexesAddition.insert(elem);
                    }
                    
                }
            }
            
            if (TempVertexesAddition.size() == 0) {
                break;
            } else {
                
                rand_vertex = rand() % TempVertexesAddition.size();
                
                vertex = TempVertexesAddition.begin();
                
                for (int i = 0; i < rand_vertex; ++i) {
                    vertex++;
                }
                
                S.Add(*vertex);
            }
            
        }
        return S;
    }
};

class Metropolis final: public VertexCoverSolver {


    public:
        Metropolis(int k, int t, bool flag);  // This is the constructor

    bool checker(const VertexCover& vertexset, Graph::Vertex& elem_to_check) const {
        bool flag = false;
        for (auto elem : vertexset) {
            if (elem == elem_to_check) {
                flag == true;
                break;
            }
        }
        
        return flag;
    }


    VertexCover Solve(int k, int t,bool flag, const Graph& graph, DebugInfo& debug_info) const {
        
        int rand_vertex = rand() % graph.AllVertices().size();
        VertexCover S(graph);
        
        auto vertex = graph.AllVertices().begin();
        
        int i = 0;
        
        for (auto elem : graph.AllVertices()){
            if (i == rand_vertex) {
                S.Add(elem);
            }
            i++;
        }
        
        while (true) {
            
            Graph::VertexSet TempVertexesAddition;
            
            for (auto vert = S.begin(); vert != S.end(); ++vert) {
                Graph::VertexSet temp = graph.AdjacentVertices(*vert);
                
                for (auto elem : temp) {
                    
                    if (checker(S, elem)) {
                        TempVertexesAddition.insert(elem);
                    }
                    
                }
            }
            
            if (TempVertexesAddition.size() == 0) {
                break;
            } else {
                
                rand_vertex = rand() % TempVertexesAddition.size();
                
                vertex = TempVertexesAddition.begin();
                
                for (int i = 0; i < rand_vertex; ++i) {
                    vertex++;
                }
                
                S.Add(*vertex);
            }
            
        }
        return S;
    }
};

Graph RandomGraph(size_t size, double edge_probability) {
    Graph graph;
    for (Graph::Vertex v = 1; v <= size; ++v) {
        graph.AddVertex(v);
    }
    for (Graph::Vertex v = 1; v <= size; ++v) {
        for (Graph::Vertex u = v + 1; u <= size; ++u) {
            if (double(rand()) / RAND_MAX <= edge_probability) {
                graph.AddEdge(v, u);
            }
        }
    }
    return graph;
}

Graph StarGraph(size_t size) {
    Graph graph;
    for (Graph::Vertex v = 2; v <= size; ++v) {
        graph.AddEdge(1, v);
    }
    return graph;
}


void TrySolver(const VertexCoverSolver& solver, const Graph& graph) {
    GraphViz(std::cout, graph);
    auto best_cost = std::numeric_limits<size_t>::max();
    size_t results = 0;
    for (int attempt = 1; attempt < 100; ++attempt) {
        DebugInfo debug_info;
        const auto vertex_cover = solver.Solve(graph, debug_info);
        auto cost = vertex_cover.Cost();
        if (cost < best_cost) {
            best_cost = cost;
            GraphViz(std::cout, vertex_cover);
            std::cout << "Trace info:\n" << debug_info << "\n";
            ++results;
        }
    }
    std::cout << "Results: " << results << std::endl;
}

int main(int argc, const char* argv[]) {
    std::cout << "Using rand seed: " << InitRandSeed(argc, argv) << "\n";
    
    const auto graph = RandomGraph(999, 0.5);
    GradientDescent gradient_descent;
    Metropolis metropolis(1, 100, true);
    TrySolver(gradient_descent, graph);
    TrySolver(metropolis, graph);
    
    return 0;
}