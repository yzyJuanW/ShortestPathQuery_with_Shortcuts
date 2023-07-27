#ifndef TREE_DECOMP_H
#define TREE_DECOMP_H

#include <iostream>
#include <sstream>
#include <map>
#include <queue>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>
#include <algorithm> 
#include "PLF.h"
using namespace std;

typedef unsigned int NodeId;
typedef unsigned int EdgeId;
typedef pair<double, pair<NodeId,NodeId>> UPshortcut; 
//typedef pair<int,NodeId> node; // (deg, vertex) -> (height, vertex)
typedef pair<NodeId, PLF> Node2PLF;
struct vEdge {
    unsigned int dst;
    PLF weights;
};


struct TreeNode {
    int height;
    int pnodeid;
    vector<NodeId> cnodeid;

    ////////  edgesPLF1: X(v) to X(v)\{v}; edgesPLF2: X(v)\{v} to X(v); 
    vector<PLF> edgesPLF1;
    vector<PLF> edgesPLF2;

    vector<NodeId> edgesId;
    vector<int> position;

    vector<NodeId> ancestors;
    // vector<pair<int,double>> CandidateUPshortcutPair;// int->size double->utility

    vector<PLF> selectShortcurts1;
    vector<PLF> selectShortcurts2;
    vector<int> shortcurtsPosition1;
    vector<int> shortcurtsPosition2;

};

vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
////////////////////////////////////////////////////////////////////////////

class TDTree
{
    public:
    int vorder = 0;
    int width = 0;
    NodeId n;
    EdgeId m;
    NodeId root;
    vector<vector<vEdge> > graph;
    vector<TreeNode> tnodes;
    vector<int> vertexOrder;
    vector<vector<NodeId> > height2nodeQ;//height->a list of vertex at the same height
    
    //for 2-hop index:
    vector<NodeId> core_vertexes; 
    vector<int> min_degre;

    vector<vector<PLF>> BuiltShortcuts1;  // All shortcuts in the tree
    vector<vector<PLF>> BuiltShortcuts2;  // All shortcuts in the tree



    /////////////////////////////////////////////////// Debug Functions ///////////////////////////////////////////////////
    void report_shortcutes(vector<map<unsigned int, PLF> > &shortcuts)
    {
        cout << "Report Shortcurts ****************************: " << endl;
        for (int ii = 0; ii < graph.size(); ii++)
        {
            auto &vv = shortcuts[ii];
            cout << "From " << ii << " to: " << endl;
            for (map<unsigned int, PLF>::iterator it = vv.begin(); it!=vv.end(); ++it)
            {
                cout << it->first <<" ,";
                cout << it->second << endl;
            }
        }      
    }

    void report_degree2nodeQ(vector<vector<NodeId> > &degree2nodeQ)
    {
        cout << " Report degree2nodeQ ****************************: " << endl;
        for (int i = 0; i < degree2nodeQ.size(); i++)
        {
            cout << "degree: " << i <<";  " << " #" << degree2nodeQ[i].size() << endl;            
        }
    }

    void report_vPosition(vector<pair<NodeId, NodeId> > &vPosition)
    {
        cout << " Report vPosition ****************************: " << endl;
        for (int i = 0; i < graph.size(); i++)
        {
            cout << "vertex:" << i << " degree: " << vPosition[i].first << " index:" << vPosition[i].second << endl;
        }
        cout << endl;
    }    

    void report_vertexorder()
    {
        cout << "Report VertexOrder: " << endl;
        for (int ii = 0; ii < graph.size(); ii++)
        {
            cout << ii << ":" << vertexOrder[ii] << ", ";
        }
        cout << endl;
    }

    ///////////////////////////////////////////////////                 ///////////////////////////////////////////////////

    void load_graph_zy()
    {
        clock_t tbegin, tend;
        tbegin = clock();

        ifstream readgraph(tgraph_path);
        assert(readgraph.is_open());
        int a,b;
        readgraph >> n >> m >> a >> b;
        cout << "vertexes: " << n << " edges: " << m <<endl;
        graph.resize(n);
        int vs, vt, weight_piece_num;
        while (readgraph >> vs >> vt >> weight_piece_num) 
        {//input edges
            vEdge ve;
            ve.dst = vt;
            for (int i = 0; i < weight_piece_num; i++)
            {
                double t,w;
                readgraph >> t >> w;
                Segment seg(t,w);
                ve.weights.f->push_back(seg);
            }

            graph[vs].push_back(ve);

        }
        
        tend = clock();
        cout << "Finish loading graph \t time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;
        // cout << "============================report graph========================" << endl;
        // for (int i = 0; i < graph.size(); i++)
        // {
        //     for (int j = 0; j < graph[i].size(); j++)
        //     {
        //         auto &e = graph[i][j];
        //         cout << i << " " << e.dst << endl;
        //         cout << e.weights << endl;
        //         cout << " -------------------------------- " << endl;
        //     }
            
        // }
        // cout << "==========================report end======================" << endl;        
    }

    void report() 
    {
        cout << "============================report========================" << endl;
        cout << "---tree index---" << endl;
        int max_height = -100;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            cout << "vid: " << i+1 << endl;
            cout << "order: " << vertexOrder[i] << endl;
            cout << "pnode: " << tnodes[i].pnodeid+1 << endl;
            cout << "height: " << tnodes[i].height << endl;
            cout << "cnodes: ";
            for (NodeId c:tnodes[i].cnodeid)cout << c+1 << ", ";
            cout << endl;
            cout << "X(node): " ;
            for (int j = 0; j < tnodes[i].edgesId.size(); j++)
            {
                cout << tnodes[i].edgesId[j] + 1 << ": " << tnodes[i].edgesPLF1[j] << "; " << tnodes[i].edgesPLF2[j]<< "; " ;
            }
            
            cout << endl;
            cout << "Position: " ;
            for (int p:tnodes[i].position)
            {
                cout << p <<", ";
            }
            cout << endl;            
            cout << "Ancestors: ";
            for(NodeId avid:tnodes[i].ancestors)
            {
                cout << avid+1 << ", ";
            }
            if(tnodes[i].height > max_height) max_height = tnodes[i].height;
            cout << endl;

            cout << "Build shortcuts: " << endl;
            for (int j = 0; j < tnodes[i].shortcurtsPosition1.size(); j++)
            {
                if (tnodes[i].shortcurtsPosition1[j] >= 0)
                {
                    cout << tnodes[i].ancestors[j]+1 << " :" << tnodes[i].selectShortcurts1[j] << endl;
                }
                
            }
            cout << endl;
            if(tnodes[i].height > max_height) max_height = tnodes[i].height;
            cout << endl;            
            cout << " --------------------------------------" << endl;
        }
        
        cout << "max height: " << max_height << endl;
        cout << "==========================report end======================" << endl;
    }    

    unsigned int report_idx_size() 
    {

        int intNum = 0, doubleNum = 0;
        intNum += vertexOrder.size();
        intNum += core_vertexes.size();

        for (int i = 0; i < graph.size(); i++)
        {
            for (int j = 0; j < graph[i].size(); j++)
            {
                auto &e = graph[i][j];
                intNum += 1;
                intNum += e.weights.f->size();
                doubleNum += 2*e.weights.f->size();
            }
            
        }

        for (int i = 0; i < tnodes.size(); i++)
        {
            intNum += 2;
            intNum += tnodes[i].cnodeid.size();
            intNum += tnodes[i].ancestors.size();

            for (int j = 0; j < tnodes[i].edgesId.size(); j++)
            {
                intNum += 2;
                intNum += tnodes[i].edgesPLF1[j].f->size();
                intNum += tnodes[i].edgesPLF2[j].f->size();
                doubleNum += tnodes[i].edgesPLF1[j].f->size()*2;
                doubleNum += tnodes[i].edgesPLF2[j].f->size()*2;
            }
            
            for (int j = 0; j<tnodes[i].selectShortcurts1.size(); j++) 
            {
                intNum += 1;
                intNum += tnodes[i].selectShortcurts1[j].f->size();
                doubleNum += tnodes[i].selectShortcurts1[j].f->size()*2;                
            }
            intNum += tnodes[i].shortcurtsPosition1.size();

            for (int j = 0; j<tnodes[i].selectShortcurts2.size(); j++) 
            {
                intNum += 1;
                intNum += tnodes[i].selectShortcurts2[j].f->size();
                doubleNum += tnodes[i].selectShortcurts2[j].f->size()*2;                
            }
            intNum += tnodes[i].shortcurtsPosition2.size();
        }
        
        unsigned int idxsize = (intNum*4 + doubleNum*8)/1000000;
        cout << "idx_size= " <<idxsize << "MB" << endl;
        return idxsize;    
    }

    void cal_parentchild_relationship_of_tree_zy() 
    {
        for (unsigned int vid = 0; vid < tnodes.size(); ++vid) 
        {
            auto &node = tnodes[vid];
            auto &edgesId = node.edgesId;

            if(edgesId.size() >0)
            {
                //cout << "For node: " << vid+1 << endl;
                NodeId neighbor = edgesId[0];
                // NodeId idx = neighbor;
                // int order = vertexOrder[neighbor];

                //NodeId idx = 0;
                int order = 1000000000;
                node.pnodeid = neighbor;

                for (int j = 0; j < edgesId.size(); j++)
                {
                    NodeId x = edgesId[j];
                    //cout << x+1 << "  " << vertexOrder[x] << endl;
                    if (vertexOrder[x] < order and vertexOrder[x] > 0)
                    {
                        //cout << x+1 << "  " << vertexOrder[x] << endl;
                        order = vertexOrder[x];
                        node.pnodeid = x;
                        //idx = x;
                        //cout << x+1 << "  " << order << endl;
                    }
                }
                //cout << "pnodeid: " << idx+1 << endl;
                //node.pnodeid = idx;
                //tnodes[idx].cnodeid.push_back(vid);
                tnodes[node.pnodeid].cnodeid.push_back(vid);
                
            }
            else
            {
                node.pnodeid = -1;
            }
        }
    }

    void calculate_height()
    {
        queue<NodeId> Q;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            if (tnodes[i].pnodeid == -1)
            {
                Q.push(i);
                root = i;
                tnodes[i].height = 1;
            }

            if(width < tnodes[i].edgesPLF1.size())
            {
                width = tnodes[i].edgesPLF1.size();
            }   
        }

        while (!Q.empty()) 
        {
            auto v = Q.front();
            Q.pop();
            for (auto a:tnodes[v].cnodeid) 
            {
                Q.push(a);
                tnodes[a].height = tnodes[v].height + 1;
            }
        }
    }

    void calculate_ancestors()
    {
        queue<NodeId> Q;
        // find tree root
        for (unsigned int vid = 0; vid < tnodes.size(); ++vid) 
        {
            if(tnodes[vid].pnodeid == -1) Q.push(vid);
        }

        while (!Q.empty()) 
        {
            NodeId v = Q.front();
            for(NodeId u:tnodes[v].cnodeid) Q.push(u);

            NodeId pvid = tnodes[v].pnodeid;
            if (pvid == -1)
            {
                Q.pop();
                continue;
            }
            for(NodeId anc:tnodes[pvid].ancestors) tnodes[v].ancestors.push_back(anc);
            tnodes[v].ancestors.push_back(pvid);

            Q.pop();
        }
        
        //cout << "Finish calculate all ancestors.. " << endl;

        // caculate positions
        for (unsigned int vid = 0; vid < tnodes.size(); ++vid) 
        {
            //cout << "Calculate Positions: " << vid <<endl;
            for (unsigned int uid: tnodes[vid].edgesId)
            {
                auto it = find(tnodes[vid].ancestors.begin(), tnodes[vid].ancestors.end(), uid);
                tnodes[vid].position.push_back(int(it - tnodes[vid].ancestors.begin()));    
            }
            
            //step: resize & initial shortcuts
            tnodes[vid].shortcurtsPosition1.resize(tnodes[vid].ancestors.size(),0);
            tnodes[vid].shortcurtsPosition2.resize(tnodes[vid].ancestors.size(),0);
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            tnodes[vid].selectShortcurts1.resize(tnodes[vid].ancestors.size(),weights);
            tnodes[vid].selectShortcurts2.resize(tnodes[vid].ancestors.size(),weights);
        }
    }

    void graph_reduction_zy()
    {
        vector<vector<NodeId> > degree2nodeQ;//degree->a list of vertex of this degree
        vector<pair<NodeId, NodeId> > vPosition(graph.size());//(degree,idx)
        for (NodeId v = 0; v < graph.size(); ++v) 
        {
            NodeId degree = graph[v].size();
            if (degree >= degree2nodeQ.size()) {
                degree2nodeQ.resize(degree + 1);
            }
            vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
            degree2nodeQ[degree].push_back(v);
        }  

        // report_degree2nodeQ(degree2nodeQ);
        // return;

        vector<map<unsigned int, PLF> > shortcuts(graph.size());
        for (NodeId s = 0; s < graph.size(); s++)
        {
            for (unsigned int i = 0; i < graph[s].size(); ++i)
            {
                auto &e = graph[s][i];
                shortcuts[s][e.dst] = e.weights;
            }
        }  

        vertexOrder.resize(graph.size(), -1);
        NodeId mindegree = 0;
        tnodes.resize(graph.size());

        vector<NodeId> vertexPreOrder = {9,11,4,1,5,10,6,8,0,7,2,3};

        for (int i = 0; i < graph.size() - 1; i++)
        {
            mindegree = 0;
            while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) 
            {
                mindegree++;
            }

            NodeId vid = degree2nodeQ[mindegree].back();
            //NodeId vid = vertexPreOrder[i];  // debug for vector<NodeId> vertexPreOrder = {9,11,4,1,5,10,6,8,0,7,2,3};
            degree2nodeQ[mindegree].pop_back();
            vertexOrder[vid] = vorder++;  

            min_degre.push_back(mindegree);
            // reduing on node vid
            cout << "  ===================================== Reducing node =====================================: " << vid+1 << endl;
            cout << "  ===================================== Mnidegree =====================================: " << mindegree << endl;
            auto &v = shortcuts[vid];

            // Reducing step 1: find nbr(vid, G_{i-1})
            vector<unsigned int> valid_neighbor_index;
            for (map<unsigned int, PLF>::iterator it = v.begin(); it!=v.end(); ++it)
            {
                // vertexOrder[id]!=-1 => vertex id has been eliminated 
                if (vertexOrder[it->first]==-1) 
                {
                    valid_neighbor_index.push_back(it->first);
                }  
            }

            // Reducing step 2: find u,w in nbr(vid, G_{i-1})
            vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), -1);
            for (unsigned int ii = 0; ii < valid_neighbor_index.size(); ++ii) 
            {
                for (unsigned int jj = ii + 1; jj < valid_neighbor_index.size(); ++jj) 
                {
                    NodeId ivid  = valid_neighbor_index[ii];
                    NodeId jvid  = valid_neighbor_index[jj];
                    
                    PLF PLFij, PLFji;
                    shortcuts[vid][jvid].compound(shortcuts[ivid][vid],PLFij,vid);
                    shortcuts[jvid][vid].compound(shortcuts[vid][ivid],PLFji,vid);
                    if (shortcuts[ivid].find(jvid) == shortcuts[ivid].end())
                    {
                        neighbor_degree_increase_cnt[ii]++;
                        neighbor_degree_increase_cnt[jj]++;
                        shortcuts[ivid][jvid] = PLFij;
                        shortcuts[jvid][ivid] = PLFji;
                    }
                    else
                    {
                        shortcuts[ivid][jvid].minimize(PLFij);
                        shortcuts[jvid][ivid].minimize(PLFji);
                    }

                }

            }

            // Reducing step 3: update vPosition of compound vertexes
            for (unsigned int i = 0; i < valid_neighbor_index.size(); ++i) 
            { 
                //if (neighbor_degree_increase_cnt[i] > 0) 
                if (neighbor_degree_increase_cnt[i] != 0) 
                {
                    //cout << "1. swap and delete " << endl;
                    int x = valid_neighbor_index[i];
                    //cout << " pair<NodeId, NodeId> &p = vPosition[x]; " << vPosition[x].first << " " << vPosition[x].second << endl;
                    pair<NodeId, NodeId> &p = vPosition[x];
                    //1. swap and delete
                    //cout << " degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back(); " << degree2nodeQ[p.first].size() << " " << p.second << endl;
                    degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
                    //cout << " vPosition[degree2nodeQ[p.first].back()].second = p.second; " << degree2nodeQ[p.first].back() << " " << p.second << endl;
                    vPosition[degree2nodeQ[p.first].back()].second = p.second;
                    //cout << " degree2nodeQ[p.first].pop_back(); " << p.first << " " << p.second << endl;
                    degree2nodeQ[p.first].pop_back();
                    //2. place in a new position
                    //cout << "place in a new position" << endl;
                    //p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
                    //cout << "2. place in a new position " << endl;
                    //cout << " test :  p.first += neighbor_degree_increase_cnt[i];"<< p.first << " " << neighbor_degree_increase_cnt[i] << endl;
                    p.first += neighbor_degree_increase_cnt[i];
                    if (p.first == 0)
                    {
                        continue;
                    }
                    
                    //cout << " after :  p.first += neighbor_degree_increase_cnt[i];"<< p.first << endl;
                    if (p.first >= degree2nodeQ.size()) 
                    {
                        degree2nodeQ.resize(p.first + 1);
                    }
                    //cout << "mindegree = min(mindegree, p.first);" << endl;
                    //mindegree = min(mindegree, p.first);
                    //cout << "p.second = degree2nodeQ[p.first].size();" << endl;
                    //cout << "p.first: " << p.first << " " << "mindegree: " << mindegree << endl;
                    //cout << "degree2nodeQ.size(): " << degree2nodeQ.size() << endl;
                    p.second = degree2nodeQ[p.first].size();
                    //cout << " degree2nodeQ[p.first].push_back(x);" << endl;
                    degree2nodeQ[p.first].push_back(x);
                }
                //cout << " ------------------ finish neighbor: ------------ " << i << " " << valid_neighbor_index.size() << endl;
            }

            // Reducing step 4: building tnodes

            // cout << "Reducing step 4: building tnodes" << endl;
            // cout << valid_neighbor_index.size() << endl;
            for(auto i:valid_neighbor_index) 
            {
                // cout << i << " " << v.size() <<endl;
                // cout << v[i] << endl;
                // cout << tnodes.size() << endl;
                //tnodes[vid].edges[i] = v[i];
                // tnodes[vid].edgesId.push_back(i);
                // tnodes[vid].edges.push_back(v[i]);
                tnodes[vid].edgesId.push_back(i);
                tnodes[vid].edgesPLF1.push_back(v[i]);
                tnodes[vid].edgesPLF2.push_back(shortcuts[i][vid]);

                //cout << "Reduction on node :" << vid+1 << " ->" << i+1 << endl;
                // cout << tnodes[vid].edges[i] << " &&&&&&&&&&&&&&&& " << endl;
                //tnodes[vid].edges.insert(std::move(v.idx(i)));
            }
            if (tnodes[vid].edgesPLF1.size()==0) {
                //cout << "if (tnodes[vid].edges.size()==0)" << endl; 
                core_vertexes.emplace_back(vid);
            }          
        }
        
        //cout << "cal_parentchild_relationship_of_tree();" << endl;
        cal_parentchild_relationship_of_tree_zy();
        //cout << "Finish cal_parentchild_relationship_of_tree_zy();" << endl;
        //cout << "start calculate_height();" << endl;
        calculate_height();
        //cout << "Finish calculate_height();" << endl;
        //cout << "start calculate_ancestors();" << endl;   
        calculate_ancestors();
        //cout << "start cal_CandidateUPshortcutPair();" << endl;   
        //cal_CandidateUPshortcutPair();

    }

    void store_index()
    {
        ofstream index_write(index_path);

        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            index_write << to_string(i) << " " << to_string(vertexOrder[i]) << " " << to_string(tnodes[i].pnodeid) << " " << to_string(tnodes[i].height) << "\n";
            
            index_write << tnodes[i].cnodeid.size() << "\n";
            for (auto c:tnodes[i].cnodeid)index_write << c << "\n";

            index_write << tnodes[i].ancestors.size() << "\n";
            for (auto c:tnodes[i].ancestors)index_write << c << "\n";
           
            index_write << tnodes[i].position.size() << "\n";
            for (auto c:tnodes[i].position)index_write << c << "\n";

            index_write << tnodes[i].edgesId.size() << "\n";
            for (int j = 0; j < tnodes[i].edgesId.size(); j++)
            {                
                index_write << tnodes[i].edgesId[j] << " " << tnodes[i].edgesPLF1[j].f->size() << "\n";
                for (int k = 0; k < tnodes[i].edgesPLF1[j].f->size(); k++)
                {
                    index_write << (*(tnodes[i].edgesPLF1[j].f))[k].t << " " << (*(tnodes[i].edgesPLF1[j].f))[k].w << " " << (*(tnodes[i].edgesPLF1[j].f))[k].intv << "\n";
                }
            }

            for (int j = 0; j < tnodes[i].edgesId.size(); j++)
            {                
                index_write << tnodes[i].edgesId[j] << " " << tnodes[i].edgesPLF2[j].f->size() << "\n";
                for (int k = 0; k < tnodes[i].edgesPLF2[j].f->size(); k++)
                {
                    index_write << (*(tnodes[i].edgesPLF2[j].f))[k].t << " " << (*(tnodes[i].edgesPLF2[j].f))[k].w << " " << (*(tnodes[i].edgesPLF2[j].f))[k].intv << "\n";
                }
            }

            for (int j = 0; j < tnodes[i].selectShortcurts1.size(); j++)
            {                
                index_write << tnodes[i].shortcurtsPosition1[j] << " " << tnodes[i].selectShortcurts1[j].f->size() << "\n";
                for (int k = 0; k < tnodes[i].selectShortcurts1[j].f->size(); k++)
                {
                    index_write << (*(tnodes[i].selectShortcurts1[j].f))[k].t << " " << (*(tnodes[i].selectShortcurts1[j].f))[k].w << " " << (*(tnodes[i].selectShortcurts1[j].f))[k].intv << "\n";
                }
            }

            for (int j = 0; j < tnodes[i].selectShortcurts2.size(); j++)
            {                
                index_write << tnodes[i].shortcurtsPosition2[j] << " " << tnodes[i].selectShortcurts2[j].f->size() << "\n";
                for (int k = 0; k < tnodes[i].selectShortcurts2[j].f->size(); k++)
                {
                    index_write << (*(tnodes[i].selectShortcurts2[j].f))[k].t << " " << (*(tnodes[i].selectShortcurts2[j].f))[k].w << " " << (*(tnodes[i].selectShortcurts2[j].f))[k].intv << "\n";
                }
            }             
        }      
    }

    void restore_index()
    {
        ifstream index_read;
        index_read.open(index_path.c_str());
        if (!index_read.is_open()) {
            printf("%s does not exist\n", index_path.c_str());
            exit(0);
        }

        vertexOrder.resize(graph.size());
        tnodes.resize(graph.size());
        string line;
        string delimiter = " ";
        while(getline(index_read, line))
        {
            //cout << "vector<string> nodes = split(line, delimiter);" << endl;
            vector<string> nodes = split(line, delimiter);
            unsigned int vid = stoi(nodes[0]);
            vertexOrder[vid] = stoi(nodes[1]);
            //cout << "vertexOrder[vid] = stoi(nodes[1]);" << endl;
            tnodes[vid].pnodeid = stoi(nodes[2]);
            //cout << "tnodes[vid].pnodeid = stoi(nodes[2]);" << endl;
            tnodes[vid].height = stoi(nodes[3]);
            //cout << "tnodes[vid].height = stoi(nodes[3]);" << endl;
            //cout << "node: " << vid + 1 << " vertexOrder:" << vertexOrder[vid] << ",  " <<  "pnodeid: " << tnodes[vid].pnodeid+1 << "height: " << tnodes[vid].height << endl;

            string cnodeidNumline;
            getline(index_read, cnodeidNumline);
            //cout << "vector<string> cnodeidNumstring = split(cnodeidNumline, delimiter);" << endl;
            vector<string> cnodeidNumstring = split(cnodeidNumline, delimiter); 
            //cout << "cnodeid size :" << stoi(cnodeidNumstring[0]) << endl;
            for (int i = 0; i < stoi(cnodeidNumstring[0]); i++)
            {
                string cnodeidline;
                getline(index_read, cnodeidline);
                //cout << "vector<string> cnodeidstring = split(cnodeidline, delimiter);" << endl;
                vector<string> cnodeidstring = split(cnodeidline, delimiter);
                //cout << "cnodeid: ";
                //cout << stoi(cnodeidstring[0])+1 << endl;
                tnodes[vid].cnodeid.push_back(stoi(cnodeidstring[0]));            
            }
            //cout << endl;

            string ancestorNumline;
            getline(index_read, ancestorNumline);
            //cout << "vector<string> ancestorNumstring = split(ancestorNumline, delimiter);" << endl;
            vector<string> ancestorNumstring = split(ancestorNumline, delimiter);
            //cout << "ancester size :" << stoi(ancestorNumstring[0]) << endl;
            for (int i = 0; i < stoi(ancestorNumstring[0]); i++)
            {
                string ancestorline;
                getline(index_read, ancestorline);
                //cout << "vector<string> ancestorstring = split(ancestorline, delimiter);" << endl;
                vector<string> ancestorstring = split(ancestorline, delimiter);
                //cout << "ancestors: ";
                //cout << stoi(ancestorstring[0])+1 << endl;
                tnodes[vid].ancestors.push_back(stoi(ancestorstring[0]));
            }
            //cout << endl; 

            string positionNumline;
            getline(index_read, positionNumline);
            //cout << "vector<string> positionNumstring = split(positionNumline, delimiter);" << endl;
            vector<string> positionNumstring = split(positionNumline, delimiter);
            //cout << "X(node)/node size: " << stoi(positionNumstring[0]) << endl;
            for (int i = 0; i < stoi(positionNumstring[0]); i++)
            {
                string positionline;
                getline(index_read, positionline);
                //cout << "vector<string> positionstring = split(positionline, delimiter);" << endl;
                vector<string> positionstring = split(positionline, delimiter);
                tnodes[vid].position.push_back(stoi(positionstring[0]));
                //cout << "positions: ";
                //cout << stoi(positionstring[0])<<endl;       
            }                         
            //cout << endl;

            string edgesNumline;
            getline(index_read, edgesNumline);
            //cout << "vector<string> edgesNumlinestring = split(edgesNumline, delimiter); " << endl;
            vector<string> edgesNumlinestring = split(edgesNumline, delimiter);
            //cout << "edgesNum: " <<  stoi(edgesNumlinestring[0]) << endl;
            for (int i = 0; i < stoi(edgesNumlinestring[0]); i++)
            {
                string segsNumline;
                getline(index_read, segsNumline);
                //cout << "vector<string> segsNumlinetring = split(segsNumline, delimiter); " << endl;
                vector<string> segsNumlinetring = split(segsNumline, delimiter);
                int std =  stoi(segsNumlinetring[0]);
                tnodes[vid].edgesId.push_back(std);
                tnodes[vid].edgesPLF1.resize(i+1);
                tnodes[vid].edgesPLF2.resize(i+1);

                for (int j = 0; j < stoi(segsNumlinetring[1]); j++)
                {
                    string weightline;
                    getline(index_read, weightline);
                    //cout << "vector<string> weightsstring = split(weightline, delimiter); " << endl;
                    vector<string> weightsstring = split(weightline, delimiter);
                    double t = stod(weightsstring.at(0));
                    double w = stod(weightsstring.at(1));
                    int intv = stoi(weightsstring.at(2));
                    Segment seg(t,w,intv);
                    tnodes[vid].edgesPLF1[i].f->push_back(seg);
                }

                for (int j = 0; j < stoi(segsNumlinetring[1]); j++)
                {
                    string weightline;
                    getline(index_read, weightline);
                    //cout << "vector<string> weightsstring = split(weightline, delimiter); " << endl;
                    vector<string> weightsstring = split(weightline, delimiter);
                    double t = stod(weightsstring.at(0));
                    double w = stod(weightsstring.at(1));
                    int intv = stoi(weightsstring.at(2));
                    Segment seg(t,w,intv);
                    tnodes[vid].edgesPLF2[i].f->push_back(seg);
                }
                                
            }

            // string candidatePairNumline;
            // getline(index_read, candidatePairNumline);
            // //cout << "vector<string> edgesNumlinestring = split(edgesNumline, delimiter); " << endl;
            // vector<string> candidatePairNumlinestring = split(candidatePairNumline, delimiter);
            // //cout << "edgesNum: " <<  stoi(edgesNumlinestring[0]) << endl;
            // for (int i = 0; i < stoi(candidatePairNumlinestring[0]); i++)
            // {
            //     string pairline;
            //     getline(index_read, pairline);
            //     //cout << "vector<string> weightsstring = split(weightline, delimiter); " << endl;
            //     vector<string> pairstring = split(pairline, delimiter);
            //     tnodes[vid].CandidateUPshortcutPair.push_back(make_pair(stoi(pairstring[0]),stod(pairstring[1])));
            // }

               
        }
    }    

    pair<unsigned int, unsigned int> qLCA(NodeId s, NodeId t) 
    {
        NodeId sAncestor = s, tAncestor = t;
        while (true) 
        {
            if (tnodes[sAncestor].height > tnodes[tAncestor].height) 
            {
                if (static_cast<int>(tAncestor)==tnodes[sAncestor].pnodeid) return make_pair(tAncestor, tAncestor);
                sAncestor = static_cast<unsigned int>(tnodes[sAncestor].pnodeid);
            } 
            else if (tnodes[sAncestor].height < tnodes[tAncestor].height) 
            {
                if (static_cast<int>(sAncestor)==tnodes[tAncestor].pnodeid) return make_pair(sAncestor, sAncestor);
                tAncestor = static_cast<unsigned int>(tnodes[tAncestor].pnodeid);
            } 
            else 
            {
                if (tnodes[sAncestor].pnodeid==-1 || sAncestor==tAncestor)return make_pair(sAncestor, tAncestor);
                sAncestor = static_cast<unsigned int>(tnodes[sAncestor].pnodeid);
                tAncestor = static_cast<unsigned int>(tnodes[tAncestor].pnodeid);
            }
        }
    }

   
    void creat_index() 
    {
        graph_reduction_zy();
    }

    void ShowTreeStructure()
    {
        ofstream write_tree("tree.txt");

        int num = 117756;
        int height = 10;
        // random show num times tree, tree with height 
        for (int i = 117754; i < num; i++)
        {
            //NodeId root = rand()%(graph.size()-1);
            NodeId root = i;
            NodeId s  = root;
            queue<NodeId> que;
            cout << "------------------- Root:  -----------------" << s << endl;

            cout << " root graph edges: ";
            for (auto edges: graph[s])
            {
                cout << edges.dst << " ,";
            }
            cout << endl;

            write_tree << "------------------- Root:  -----------------" << s << "\n";
            write_tree << "----- Root.pnode:  -----" << tnodes[s].pnodeid << "\n";
            for(NodeId cnode:tnodes[s].cnodeid) que.push(cnode);
            while (!que.empty())
            {
                s  = que.front();
                cout << s  << "::  parent: " << tnodes[s].pnodeid << " X(s): ";
                write_tree << "For node " <<  s  << "->  parent: " << tnodes[s].pnodeid << " X(s): ";
                for (NodeId p:tnodes[s].edgesId)
                {
                    cout << p << ", ";
                    write_tree << p << ", ";
                } 
                cout << endl;
                write_tree << "\n";
                cout << "cnode: " ;
                write_tree << "cnode: " ;
                for (NodeId p:tnodes[s].cnodeid)
                {
                    if (tnodes[p].height <= tnodes[root].height + height)
                    {
                        cout << p <<", ";
                        write_tree << p <<", ";
                        que.push(p);
                    }
                    
                }

                //graph
                cout << endl;
                cout << "graph edges: ";
                for (auto edges: graph[s])
                {
                    cout << edges.dst << " ,";
                }
                cout << endl;

                que.pop();
                cout << endl;
                write_tree << "\n";

            }
            cout << " ---------------------------------- " << endl;
            write_tree << " ---------------------------------- " << "\n";
        }
        
    }

    ////////////////////////////////////////////////  Basic Query //////////////////////////////////////
    void travelCostFunctionLists(NodeId s, vector<vector<PLF>> &sResults1, vector<vector<PLF>> &sResults2)
    {
        for (int i = 0; i<=tnodes[s].ancestors.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[s].ancestors.size()+1, weights);
            sResults1.push_back(vec);
            sResults2.push_back(vec);
        }
        
        vector<NodeId> nodes;
        for (int v = 0; v < tnodes[s].ancestors.size(); v++)
        {
            nodes.push_back(tnodes[s].ancestors[v]);
        }
        nodes.push_back(s);

        //for (int v = 0; v < tnodes[s].ancestors.size(); v++)
        for (int v = 0; v < nodes.size(); v++)
        {
            NodeId nodev = nodes[v];

            for (int i = 0; i < tnodes[nodev].ancestors.size(); i++)
            {
                for (int j = 0; j < tnodes[nodev].edgesId.size(); j++)
                {
                    int posj = tnodes[nodev].position[j];

                    if (posj > i)
                    {
                        PLF PLFvji;
                        sResults1[posj][i].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }
                        

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults2[posj][i],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    else
                    {
                        PLF PLFvji;
                        sResults2[i][posj].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults1[i][posj],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    
                }
            }
        }
    }

    void basicQueryFunction(NodeId s, NodeId d, double t, PLF &costFunction)
    {
        vector<vector<PLF>> sResults1, sResults2, dResults1, dResults2; 
        travelCostFunctionLists(s, sResults1, sResults2);
        travelCostFunctionLists(d, dResults1, dResults2);

        NodeId lca = qLCA(s,d).first;

        auto itsd = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), d);
        auto itds = find(tnodes[d].ancestors.begin(), tnodes[d].ancestors.end(), s);
        if(itsd != tnodes[s].ancestors.end()) // d is an ancestor of s
        {
            int indexsd = itsd - tnodes[s].ancestors.begin();
            costFunction = sResults1[sResults1.size()-1][indexsd];
            
        }
        else if (itds != tnodes[d].ancestors.end())  // s is an ancestor of d
        {
            int indexds = itds - tnodes[d].ancestors.begin();
            costFunction = dResults2[dResults2.size()-1][indexds];
        }
        else
        {
            auto it = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), lca);
            int indexlca = it - tnodes[s].ancestors.begin();
            dResults2[dResults2.size()-1][indexlca].compound(sResults1[sResults1.size()-1][indexlca],costFunction,lca);
            for (int j:tnodes[lca].position)
            {
                PLF PLFlca;
                dResults2[dResults2.size()-1][j].compound(sResults1[sResults1.size()-1][j],PLFlca,lca);
                costFunction.minimize(PLFlca);
            }   
        }
                    
    }

    void Subgraph(NodeId s, unordered_map<NodeId, unordered_map<NodeId, PLF>> &subgraph)
    {
        for (NodeId anc:tnodes[s].ancestors)
        {
            for (int i = 0; i < tnodes[anc].edgesId.size(); i++)
            {
                subgraph[anc][tnodes[anc].edgesId[i]] = tnodes[anc].edgesPLF1[i];
                subgraph[tnodes[anc].edgesId[i]][anc] = tnodes[anc].edgesPLF2[i];
            }
            
        }
        
        for (int i = 0; i < tnodes[s].edgesId.size(); i++)
        {
            subgraph[s][tnodes[s].edgesId[i]] = tnodes[s].edgesPLF1[i];
            subgraph[tnodes[s].edgesId[i]][s] = tnodes[s].edgesPLF2[i];
        }
        
    }

    double basicQueryCost(NodeId s, NodeId d, double t)
    {
        //cout << " basicQueryCost " << endl;
        unordered_map<NodeId, unordered_map<NodeId, PLF>> sSubgraph, dSubgraph;

        Subgraph(s, sSubgraph);
        Subgraph(d, dSubgraph);

        unordered_map<NodeId, unordered_map<NodeId, PLF>> sdGraph;
        for (unordered_map<NodeId, unordered_map<NodeId, PLF>>::iterator u = sSubgraph.begin(); u!=sSubgraph.end(); u++)
        {
            for (unordered_map<NodeId, PLF>::iterator v = sSubgraph[u->first].begin(); v!=sSubgraph[u->first].end(); v++)
            {
                sdGraph[u->first][v->first] = sSubgraph[u->first][v->first];
            }
        }
        for (unordered_map<NodeId, unordered_map<NodeId, PLF>>::iterator u = dSubgraph.begin(); u!=dSubgraph.end(); u++)
        {
            for (unordered_map<NodeId, PLF>::iterator v = dSubgraph[u->first].begin(); v!=dSubgraph[u->first].end(); v++)
            {
                sdGraph[u->first][v->first] = dSubgraph[u->first][v->first];
            }
        }

        priority_queue<pair<double, NodeId>, vector<pair<double, NodeId>>, greater<pair<double, NodeId>> > Q;
        unordered_map<NodeId, double> dist;

        Q.push(make_pair(t, s));
        dist[s] = t;
        while (!Q.empty())
        {
            NodeId u  = Q.top().second;
            // double tim = QS.top().first;
            Q.pop();

            for (unordered_map<NodeId, PLF>::iterator it = sdGraph[u].begin(); it!=sdGraph[u].end(); it++)
            {
                //auto checked = find(distS.begin(), distS.end(), it->first);
                auto checked = dist.find(it->first);
                if (checked == dist.end())
                {
                    dist[it->first] = dist[u] + it->second.dpt2arr(dist[u]);
                    Q.push(make_pair(dist[it->first], it->first));
                }
                else if (dist[it->first] > dist[u] + it->second.dpt2arr(dist[u]))
                {
                    dist[it->first] = dist[u] + it->second.dpt2arr(dist[u]);
                    //Q.push(make_pair(it->first, dist[it->first]));
                    Q.push(make_pair(dist[it->first], it->first));
                } 
            }
        }
        return dist[d];
    }

    void testQuery(NodeId s)
    {
        vector<vector<PLF>> sResults1; 
        vector<vector<PLF>> sResults2; 

        travelCostFunctionLists(s, sResults1, sResults2);

        // for (int i = 0; i < tnodes[s].ancestors.size(); i++)
        // {
        //     for (int j = 0; j < i; j++)
        //     {
        //         cout << tnodes[s].ancestors[i] + 1 << " -> " << tnodes[s].ancestors[j] + 1 << endl;
        //         cout << sResults1[i][j] << endl;
        //         cout << sResults2[i][j] << endl;
        //         cout << " ------------------------------ " << endl;
        //     }
            
        //     cout << " ******************************************** " << tnodes[s].ancestors[i]+1 << " ******************************************** "<<endl;
        // }
        
    }

    ////////////////////////////////////////////////  Select shortcuts //////////////////////////////////////
    void shortcutsFunctionLists(NodeId s, vector<vector<PLF>> &sResults1, vector<vector<PLF>> &sResults2)
    {
        for (int i = 0; i<=tnodes[s].ancestors.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[s].ancestors.size()+1, weights);
            sResults1.push_back(vec);
            sResults2.push_back(vec);
        }
        
        vector<NodeId> nodes;
        for (int v = 0; v < tnodes[s].ancestors.size(); v++)
        {
            nodes.push_back(tnodes[s].ancestors[v]);
        }
        nodes.push_back(s);

        for (int v = 0; v < nodes.size(); v++)
        {
            NodeId nodev = nodes[v];

            for (int i = 0; i < tnodes[nodev].ancestors.size(); i++)
            {
                if (BuiltShortcuts1[nodev][i].dpt2wgt(0) > 0)
                {
                    sResults1[v][i] = BuiltShortcuts1[nodev][i];
                    sResults2[v][i] = BuiltShortcuts2[nodev][i];
                    continue;
                }
                

                for (int j = 0; j < tnodes[nodev].edgesId.size(); j++)
                {
                    int posj = tnodes[nodev].position[j];

                    if (posj > i)
                    {
                        PLF PLFvji;
                        sResults1[posj][i].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }
                        

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults2[posj][i],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    else
                    {
                        PLF PLFvji;
                        sResults2[i][posj].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults1[i][posj],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    
                }
            }
        }
    }    
    
    void buildShortcutsSet()
    {
        int shortcutNumLimit = 0; 

        // initialize all shortcuts 
        for (int i = 0; i<tnodes.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[i].ancestors.size(), weights);
            BuiltShortcuts1.push_back(vec);
            BuiltShortcuts2.push_back(vec);
        }
        
        vector<NodeId> LeafNodes;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            if(tnodes[i].cnodeid.size() == 0) {LeafNodes.push_back(i);}
        }
        //cout << "Leafnodes size: " << LeafNodes.size() << endl;

        for (NodeId v: LeafNodes)
        {   
            cout << " shortcutNumLimit " << shortcutNumLimit <<endl;
            vector<vector<PLF>> sResults1, sResults2; 
            shortcutsFunctionLists(v, sResults1, sResults2);
            //cout << "Finish building shortests for :" << v+1 << endl;

            for (int i = 0; i < tnodes[v].ancestors.size(); i++)
            {
                NodeId nodei = tnodes[v].ancestors[i];
                for (int j = 0; j < tnodes[nodei].ancestors.size(); j++)
                {
                    if (BuiltShortcuts1[nodei][j].dpt2wgt(0) == 0)
                    {
                        shortcutNumLimit+=2;
                        if (shortcutNumLimit>=MAXShortcutNUM){break;}
                        BuiltShortcuts1[nodei][j] = sResults1[i][j];
                        BuiltShortcuts2[nodei][j] = sResults2[i][j];
                    }
                }
                
            }
            BuiltShortcuts1[v] = sResults1[sResults1.size()-1];
            BuiltShortcuts2[v] = sResults2[sResults2.size()-1];
            if (shortcutNumLimit>=MAXShortcutNUM){break;}
        }

    }

    ///////////////////////////////////// selection algorithm 1: dp
    //vector<vector<int>> dpSelection(int N) // NodeId : ancestoer1ID...
    void dpSelection(int N) // NodeId : ancestoer1ID...
    {
        vector<pair<NodeId, int>> Candidate;
        vector<int> weights;
        vector<double> values;

        cout << "dpSelection: " << endl;

        for (unsigned int vid = 0; vid < tnodes.size(); ++vid)
        {
            auto &node = tnodes[vid];
            for (int i = 0; i < node.ancestors.size() ; i++)
            {
                Candidate.push_back(make_pair(vid,i));    //BuiltShortcuts1[][]
                // int weight = int((BuiltShortcuts1[vid][i].f->size()+BuiltShortcuts2[vid][i].f->size()));
                // weights.push_back(weight);
                weights.push_back(1);

                double utility = 0;
                double lca_size = 0, h = tnodes[vid].height - tnodes[tnodes[vid].ancestors[i]].height;
                queue<NodeId> LCAQ;
                for (NodeId lcacnode: tnodes[tnodes[vid].ancestors[i]].cnodeid)
                {
                    //if(find(tnodes[tnodes[vid].ancestors[i]].cnodeid.begin(), tnodes[tnodes[vid].ancestors[i]].cnodeid.end(), lcacnode) == tnodes[tnodes[vid].ancestors[i]].cnodeid.end()) LCAQ.push(lcacnode);
                    if(find(tnodes[vid].ancestors.begin(), tnodes[vid].ancestors.end(), lcacnode) == tnodes[vid].ancestors.end()) LCAQ.push(lcacnode);
                }
                while (!LCAQ.empty())
                {
                    NodeId v = LCAQ.front();
                    for(NodeId u:tnodes[v].cnodeid) LCAQ.push(u);                    
                    lca_size++;
                    LCAQ.pop();
                }
                lca_size++;
                utility = h*width*(lca_size/n*10000);

                values.push_back(utility);
            
            }   
        }
        cout << "Candidate.size(): " <<Candidate.size() << endl;


        int num = values.size();
        //cout << "1111111" << endl;
        vector<vector<int>> dp;
        for (int i = 0; i <= num; i++)
        {
            vector<int> vec;
            vec.resize(N+1,0);
            dp.push_back(vec);
        }
        //int dp[num+1][N+1];
        //cout << "1111111" << endl;
        for (int i = 0; i <= num; i++)
        {
            //cout << "i: " << i << endl;
            for (int w = 0; w <= N; w++)
            {
                //cout << "w: " << w<< endl;
                if (i == 0 || w==0)
                {
                    //cout << "if (i == 0 || w==0)" << endl;
                    dp[i][w] = 0;
                }
                else if(weights[i - 1] <= w)
                {
                    //cout << "else if(weights[i - 1] <= w)" << endl;
                    if(values[i - 1] + dp[i - 1][w - weights[i - 1]] > dp[i - 1][w]) {dp[i][w] =values[i - 1] + dp[i - 1][w - weights[i - 1]];}
                    else {dp[i][w] =dp[i - 1][w]; }

                }
                else 
                {
                    //cout << "else " << endl;
                    dp[i][w] =dp[i - 1][w];
                }
            }
            
        }
        
        vector<pair<NodeId, int>> selection;
        int res = dp[num][N];
        int w = N;
        for (int i = n; i>0 && res>0; i--)
        {
            if (res == dp[i-1][w])
            {
                continue;
            }
            else
            {
                selection.push_back(make_pair(Candidate[i-1].first, Candidate[i-1].second));
                res = res - values[i-1];
                w = w - weights[i-1];
            }
            
        }

        for (int i = 0; i < selection.size(); i++)
        {
            NodeId node = selection[i].first;
            int anc = selection[i].second;

            tnodes[node].shortcurtsPosition1[anc] = 1;
            tnodes[node].selectShortcurts1[anc] = BuiltShortcuts1[node][anc];

            tnodes[node].shortcurtsPosition2[anc] = 1;
            tnodes[node].selectShortcurts2[anc] = BuiltShortcuts2[node][anc];
        }
        cout << "Selection shortcut size: " << selection.size() << endl;
    }
    
    //////////////////////////////////// selection algorithm 2:  approximation algorithm
    void approSelection(int N)
    {
        cout << "start approSelection : " << endl;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            int height = tnodes[i].height;
            if (height >= height2nodeQ.size()) {
                height2nodeQ.resize(height + 1);
            }
            height2nodeQ[height].push_back(i);            
        }

        vector<pair<NodeId, int>> selection;
        int weights = 0;
        while (!height2nodeQ.empty())        
        {
            for (vector<NodeId>::iterator a = height2nodeQ.back().begin(); a!=height2nodeQ.back().end(); a++)
            {
                NodeId vid = (*a);
                auto &node = tnodes[vid];
                //cout << "check tree node :" << vid+1 <<endl;
                for (int j = 0; j < node.ancestors.size(); j++)
                {
                    if(BuiltShortcuts1[vid][j].dpt2arr(0) > 0)
                    {
                        // int weight = int((BuiltShortcuts1[vid][j].f->size()+BuiltShortcuts2[vid][j].f->size())/n);
                        // if(weight <= 0) weight = 1;
                        int weight = 1;
                        selection.push_back(make_pair(vid,j));

                        tnodes[vid].shortcurtsPosition1[j] = 1;
                        tnodes[vid].selectShortcurts1[j] = BuiltShortcuts1[vid][j];

                        tnodes[vid].shortcurtsPosition2[j] = 1;
                        tnodes[vid].selectShortcurts2[j] = BuiltShortcuts2[vid][j];  
                    

                        weights += weight;
                    }
                    if(weights >= N) break;
                }

                BuiltShortcuts1[vid].clear();
                BuiltShortcuts2[vid].clear();

                if(weights >= N) break;                
            }

            height2nodeQ.pop_back();
        }

        // for (int i = 0; i < selection.size(); i++)
        // {
        //     NodeId node = selection[i].first;
        //     int anc = selection[i].second;

        //     tnodes[node].shortcurtsPosition1[anc] = 1;
        //     tnodes[node].selectShortcurts1[anc] = BuiltShortcuts1[node][anc];

        //     tnodes[node].shortcurtsPosition2[anc] = 1;
        //     tnodes[node].selectShortcurts2[anc] = BuiltShortcuts2[node][anc];
        // }
        BuiltShortcuts1.clear();
        BuiltShortcuts2.clear();      
        cout << "Selection shortcut size: " << selection.size() << endl;
    }   

    void travelcostFunctionListsBasedShortcuts(NodeId s, vector<vector<PLF>> &sResults1, vector<vector<PLF>> &sResults2)
    {
        for (int i = 0; i<=tnodes[s].ancestors.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[s].ancestors.size()+1, weights);
            sResults1.push_back(vec);
            sResults2.push_back(vec);
        }
        
        vector<NodeId> nodes;
        for (int v = 0; v < tnodes[s].ancestors.size(); v++)
        {
            nodes.push_back(tnodes[s].ancestors[v]);
        }
        nodes.push_back(s);

        for (int v = 0; v < nodes.size(); v++)
        {
            NodeId nodev = nodes[v];

            for (int i = 0; i < tnodes[nodev].ancestors.size(); i++)
            {
                if (tnodes[nodev].shortcurtsPosition1[i] == 1 )
                {
                    sResults1[v][i] = tnodes[nodev].selectShortcurts1[i];
                    sResults2[v][i] = tnodes[nodev].selectShortcurts2[i];
                    continue;
                }
                
                for (int j = 0; j < tnodes[nodev].edgesId.size(); j++)
                {
                    int posj = tnodes[nodev].position[j];

                    if (posj > i)
                    {
                        PLF PLFvji;
                        sResults1[posj][i].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }
                        

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults2[posj][i],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    else
                    {
                        PLF PLFvji;
                        sResults2[i][posj].compound(tnodes[nodev].edgesPLF1[j], PLFvji, tnodes[nodev].edgesId[j]);
                        if (sResults1[v][i].dpt2wgt(0) == 0)
                        {
                            sResults1[v][i] = PLFvji;
                        }
                        else
                        {
                            sResults1[v][i].minimize(PLFvji);
                        }

                        PLF PLFijv;
                        tnodes[nodev].edgesPLF2[j].compound(sResults1[i][posj],PLFijv,tnodes[nodev].edgesId[j]);
                        if (sResults2[v][i].dpt2wgt(0) == 0)
                        {
                            sResults2[v][i] = PLFijv;
                        }
                        else
                        {
                            sResults2[v][i].minimize(PLFijv);
                        }
                    }
                    
                }
            }
        }
    } 

    void shortcutQueryFuncion(NodeId s, NodeId d, double t, PLF &costFunction)
    {
        vector<vector<PLF>> sResults1, sResults2, dResults1, dResults2; 
        NodeId lca = qLCA(s,d).first;

        auto itsd = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), d);
        auto itds = find(tnodes[d].ancestors.begin(), tnodes[d].ancestors.end(), s);
        if(itsd != tnodes[s].ancestors.end()) // d is an ancestor of s
        {
            int indexsd = itsd - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[indexsd] == 1)
            {
                //cout << "SHORTCUTS" << endl;
                costFunction = tnodes[s].selectShortcurts1[indexsd];
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                travelcostFunctionListsBasedShortcuts(s, sResults1, sResults2);
                costFunction = sResults1[sResults1.size()-1][indexsd];
            }
            
        }
        else if (itds != tnodes[d].ancestors.end())  // s is an ancestor of d
        {
            int indexds = itds - tnodes[d].ancestors.begin();

            if (tnodes[d].shortcurtsPosition2[indexds] == 1)
            {
                //cout << "SHORTCUTS" << endl;
                costFunction = tnodes[d].selectShortcurts2[indexds];
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                travelcostFunctionListsBasedShortcuts(d, dResults1, dResults2);
                costFunction = dResults2[dResults2.size()-1][indexds];
            }
        }
        else
        {
            int shorctcutFLAG = 1;
            auto it = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), lca);
            int lcaPosition = it - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[lcaPosition] == 0 or tnodes[d].shortcurtsPosition2[lcaPosition] == 0) {shorctcutFLAG =0;}
            for (int j: tnodes[lca].position)
            {
                if(shorctcutFLAG ==0){break;}
                if (tnodes[s].shortcurtsPosition1[j] == 0 or tnodes[d].shortcurtsPosition2[j] == 0)
                {
                    shorctcutFLAG = 0;
                    break;
                }
            }


            if (shorctcutFLAG == 1)
            {
                //cout << "SHORTCUTS" << endl;
                tnodes[d].selectShortcurts2[lcaPosition].compound(tnodes[s].selectShortcurts1[lcaPosition],costFunction,lca);
                for (int j: tnodes[lca].position)
                {
                    PLF PLFsd;
                    tnodes[d].selectShortcurts2[j].compound(tnodes[s].selectShortcurts1[j],PLFsd,tnodes[lca].ancestors[j]);
                    costFunction.minimize(PLFsd);                    
                }
                
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                travelcostFunctionListsBasedShortcuts(s, sResults1, sResults2);
                travelcostFunctionListsBasedShortcuts(d, dResults1, dResults2);

                dResults2[dResults2.size()-1][lcaPosition].compound(sResults1[sResults1.size()-1][lcaPosition],costFunction,lca);
                for (int j:tnodes[lca].position)
                {
                    PLF PLFlca;
                    dResults2[dResults2.size()-1][j].compound(sResults1[sResults1.size()-1][j],PLFlca,tnodes[lca].ancestors[j]);
                    costFunction.minimize(PLFlca);
                }          
            }
  
        }
    }

    double shortcutQueryCost(NodeId s, NodeId d, double t)
    {
        double cost = DE_W;

        NodeId lca = qLCA(s,d).first;

        auto itsd = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), d);
        auto itds = find(tnodes[d].ancestors.begin(), tnodes[d].ancestors.end(), s);
        if(itsd != tnodes[s].ancestors.end()) // d is an ancestor of s
        {
            //cout << 1111111111111 << endl;
            int indexsd = itsd - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[indexsd] == 1){cost = tnodes[s].selectShortcurts1[indexsd].dpt2arr(t); return cost;}
            else{cost = basicQueryCost(s,d,t); return cost;}
        }
        if(itds != tnodes[d].ancestors.end()) // s is an ancestor of d
        {
            //cout << 222222222222222 << endl;
            int indexds = itds - tnodes[d].ancestors.begin();
            //cout << indexds << endl;
            if (tnodes[d].shortcurtsPosition2[indexds] == 1){cost = tnodes[d].selectShortcurts2[indexds].dpt2arr(t); return cost;}
            else{cost = basicQueryCost(s,d,t); return cost;}
        }
        else
        {
            //cout << 3333333333333333333 << endl;
            int shorctcutFLAG = 1;
            auto it = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), lca);
            int lcaPosition = it - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[lcaPosition] == 1 and tnodes[d].shortcurtsPosition2[lcaPosition] == 1)
            {
                double cost1 = t + tnodes[s].selectShortcurts1[lcaPosition].dpt2arr(t);
                cost  = cost1 + tnodes[d].selectShortcurts2[lcaPosition].dpt2arr(cost1);
            }
            else{shorctcutFLAG =0;}
            for (int j: tnodes[lca].position)
            {
                if (tnodes[s].shortcurtsPosition1[j] == 1 and tnodes[d].shortcurtsPosition2[j] == 1)
                {
                    double cost1 = t + tnodes[s].selectShortcurts1[lcaPosition].dpt2arr(t);
                    double cost2 = cost1 + tnodes[d].selectShortcurts2[lcaPosition].dpt2arr(cost1);
                    if(cost2<cost){cost = cost2;}
                }
                else{shorctcutFLAG =0;}
            }

            if (shorctcutFLAG == 1) // all shortcuts are built
            {
                //cout << "SHORTCUTS" << endl;
                return cost;
            }
            else   //  "cost" is the upperbound to implement the Dijkstra
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                unordered_map<NodeId, unordered_map<NodeId, PLF>> sSubgraph, dSubgraph;

                Subgraph(s, sSubgraph);
                Subgraph(d, dSubgraph);

                unordered_map<NodeId, unordered_map<NodeId, PLF>> sdGraph;
                for (unordered_map<NodeId, unordered_map<NodeId, PLF>>::iterator u = sSubgraph.begin(); u!=sSubgraph.end(); u++)
                {
                    for (unordered_map<NodeId, PLF>::iterator v = sSubgraph[u->first].begin(); v!=sSubgraph[u->first].end(); v++)
                    {
                        sdGraph[u->first][v->first] = sSubgraph[u->first][v->first];
                    }
                }
                for (unordered_map<NodeId, unordered_map<NodeId, PLF>>::iterator u = dSubgraph.begin(); u!=dSubgraph.end(); u++)
                {
                    for (unordered_map<NodeId, PLF>::iterator v = dSubgraph[u->first].begin(); v!=dSubgraph[u->first].end(); v++)
                    {
                        sdGraph[u->first][v->first] = dSubgraph[u->first][v->first];
                    }
                }

                priority_queue<pair<double, NodeId>, vector<pair<double, NodeId>>, greater<pair<double, NodeId>> > Q;
                unordered_map<NodeId, double> dist;

                Q.push(make_pair(t, s));
                dist[s] = t;
                while (!Q.empty())
                {
                    NodeId u  = Q.top().second;
                    Q.pop();

                    for (unordered_map<NodeId, PLF>::iterator it = sdGraph[u].begin(); it!=sdGraph[u].end(); it++)
                    {
                        auto checked = dist.find(it->first);
                        if (checked == dist.end())
                        {
                            dist[it->first] = dist[u] + it->second.dpt2arr(dist[u]);
                            if (dist[it->first]<cost)
                            {
                                Q.push(make_pair(dist[it->first], it->first));
                            }
                            
                        }
                        else if (dist[it->first] > dist[u] + it->second.dpt2arr(dist[u]))
                        {
                            dist[it->first] = dist[u] + it->second.dpt2arr(dist[u]);
                            Q.push(make_pair(dist[it->first], it->first));
                        } 
                    }
                }
                return dist[d];
            }         
        }
    }

    ////////////////////////////////////////////////// Build shortcuts and query based on top->bottom manner
    void buildShortcutSetUp2Bottom()
    {
        int shorctnum = 0;
        // initialize all shortcuts 
        for (int i = 0; i<tnodes.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[i].ancestors.size(), weights);
            BuiltShortcuts1.push_back(vec);
            BuiltShortcuts2.push_back(vec);
            //cout << "Node: "  << i+1 << " ancestor size: " << tnodes[i].ancestors.size() << endl;
        }
        cout << "Root: " << root +1 <<endl;
        auto &rootNode = tnodes[root];
        queue<NodeId> Q;
        for (auto nodeid:rootNode.cnodeid)
        {
            BuiltShortcuts1[nodeid][0] = tnodes[nodeid].edgesPLF1[0];
            BuiltShortcuts2[nodeid][0] = tnodes[nodeid].edgesPLF2[0];
            shorctnum +=2;
            for(auto cnode:tnodes[nodeid].cnodeid){Q.push(cnode);}
        }
        
        while (!Q.empty())
        {
            NodeId v = Q.front();
            Q.pop();
            //cout << "Build shortcuts for:" << v+1 <<"-------------------------------------------------------- " << endl;
            for (int i = 0; i < tnodes[v].ancestors.size(); i++)
            {
                //cout << " -> " << i << "-th ancestor:"<< tnodes[v].ancestors[i]+1  << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<< endl;
                //PLF PLFvi, PLFiv;
                Segment seg(0,0);
                PLF PLFvi, PLFiv;
                PLFvi.f->push_back(seg);
                PLFiv.f->push_back(seg);

                for (int j = 0; j< tnodes[v].edgesId.size(); j++)
                {
                    //cout << tnodes[v].edgesId[j]+1 << "  " << tnodes[v].ancestors[i]+1 << endl;
                    //cout << "through: " << tnodes[v].edgesId[j]+1 << endl;
                    PLF PLFvj, PLFji, PLFvji; // v->j->i
                    PLF PLFij, PLFjv, PLFijv; // i->j->v
                    
                    if(tnodes[v].position[j] == i)
                    {
                        //cout << "case 3: j == i *****************************" << endl; 
                        PLFvji = tnodes[v].edgesPLF1[j];
                        PLFijv = tnodes[v].edgesPLF2[j];

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            // cout << "case 3.1 : i-th shortcut has not been built" << endl;
                            // cout << "shortcut PLFvi: " << PLFvji << endl;
                            // cout << "shortcut PLFiv: " << PLFijv << endl;
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;                            
                        }
                        else
                        {
                            // cout << "case 3.2 : i-th shortcut has been built" << endl;
                            // cout << "shortcut PLFvi before : " << PLFvji << endl;
                            // cout << "shortcut PLFiv before : " << PLFijv << endl;
                            //cout << "PLFvi.minimize(PLFvji);" << endl;
                            PLFvi.minimize(PLFvji);
                            //cout << "PLFiv.minimize(PLFijv);" << endl;
                            PLFiv.minimize(PLFijv);
                            // cout << "new shortcut PLFvi: " << PLFvji << endl;
                            // cout << "new shortcut PLFiv: " << PLFijv << endl;                              
                            // cout << "shortcut PLFvi after : " << PLFvi << endl;
                            // cout << "shortcut PLFiv after : " << PLFiv << endl;                               
                        }
                    }
                    else
                    {
                        if (tnodes[v].position[j] > i)
                        {
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = BuiltShortcuts1[tnodes[v].edgesId[j]][i];

                            PLFij = BuiltShortcuts2[tnodes[v].edgesId[j]][i];
                            PLFjv = tnodes[v].edgesPLF2[j];
                            
                            // cout << "case 1: j>i: v->j->i *****************************" << endl;
                            // cout << "PLFvj : " << PLFvj << endl; 
                            // cout << "PLFji : " << PLFji << endl; 
                            // cout << "PLFij : " << PLFij << endl; 
                            // cout << "PLFjv : " << PLFjv << endl; 
                        }
                        else
                        {
                            //cout << "tnodes[v].position[j] < i" << endl;  
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = BuiltShortcuts2[tnodes[v].ancestors[i]][tnodes[v].position[j]];

                            PLFij = BuiltShortcuts1[tnodes[v].ancestors[i]][tnodes[v].position[j]];
                            PLFjv = tnodes[v].edgesPLF2[j];

                            // cout << "case 2: j < i: v->j->i *****************************" << endl;
                            // cout << "PLFvj : " << PLFvj << endl; 
                            // cout << "PLFji : " << PLFji << endl; 
                            // cout << "PLFij : " << PLFij << endl; 
                            // cout << "PLFjv : " << PLFjv << endl;              
                        }
                        PLFji.compound(PLFvj,PLFvji,tnodes[v].edgesId[j]);
                        PLFjv.compound(PLFij,PLFijv,tnodes[v].edgesId[j]);

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;
                            continue;
                        }
                        PLFvi.minimize(PLFvji);
                        //cout << "PLFiv.minimize(PLFijv);" << endl;
                        PLFiv.minimize(PLFijv); 
                    }   
                    

                }
                // cout << PLFvi << endl;
                // cout << PLFiv << endl;
                BuiltShortcuts1[v][i] = PLFvi;
                BuiltShortcuts2[v][i] = PLFiv;
                shorctnum+=2; 
            }
            for(auto cnode:tnodes[v].cnodeid){Q.push(cnode);}
        }
        cout << "build shortcut: " << shorctnum << endl;
    }

    void shortcutsForNode(NodeId x, vector<vector<PLF>> &shortcuts1, vector<vector<PLF>> &shortcuts2) // shortcuts between v and all its ancestores
    {
        for (int i = 0; i<=tnodes[x].ancestors.size(); i++)
        {
            Segment seg(0,0);
            PLF weights;
            weights.f->push_back(seg);
            vector<PLF> vec;
            vec.resize(tnodes[x].ancestors.size()+1, weights);
            shortcuts1.push_back(vec);
            shortcuts2.push_back(vec);
        }

        for (int i = 0; i< tnodes[x].ancestors.size(); i++)
        {
            NodeId ancx = tnodes[x].ancestors[i];
            for (int j = 0; j < tnodes[ancx].ancestors.size(); j++)
            {
                //if (tnodes[ancx].shortcurtsPosition1[j] == 1)
                //{
                shortcuts1[i][j] = tnodes[ancx].selectShortcurts1[j];
                shortcuts2[i][j] = tnodes[ancx].selectShortcurts2[j];
                //}   
            }    
        }

        vector<NodeId> nodes;
        for (int v = 0; v < tnodes[x].ancestors.size(); v++)
        {
            nodes.push_back(tnodes[x].ancestors[v]);
        }
        nodes.push_back(x);

        for (int k = 0; k < nodes.size(); k++)
        {
            NodeId v = nodes[k];

            for (int i = 0; i < tnodes[v].ancestors.size(); i++)
            {

                if (shortcuts1[k][i].dpt2arr(0) > 0){continue;}
                
                Segment seg(0,0);
                PLF PLFvi, PLFiv;
                PLFvi.f->push_back(seg);
                PLFiv.f->push_back(seg);

                for (int j = 0; j< tnodes[v].edgesId.size(); j++)
                {

                    PLF PLFvj, PLFji, PLFvji; // v->j->i
                    PLF PLFij, PLFjv, PLFijv; // i->j->v
                    
                    if(tnodes[v].position[j] == i)
                    {
                        PLFvji = tnodes[v].edgesPLF1[j];
                        PLFijv = tnodes[v].edgesPLF2[j];

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;                            
                        }
                        else
                        {
                            PLFvi.minimize(PLFvji);
                            PLFiv.minimize(PLFijv);
                        }
                    }
                    else
                    {
                        if (tnodes[v].position[j] > i)
                        {
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = shortcuts1[tnodes[v].position[j]][i];

                            PLFij = shortcuts2[tnodes[v].position[j]][i];
                            PLFjv = tnodes[v].edgesPLF2[j];
                        
                        }
                        else
                        {
                            //cout << "tnodes[v].position[j] < i" << endl;  
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = shortcuts2[i][tnodes[v].position[j]];

                            PLFij = shortcuts1[i][tnodes[v].position[j]];
                            PLFjv = tnodes[v].edgesPLF2[j];
            
                        }
                        PLFji.compound(PLFvj,PLFvji,tnodes[v].edgesId[j]);
                        PLFjv.compound(PLFij,PLFijv,tnodes[v].edgesId[j]);

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;
                            continue;
                        }
                        PLFvi.minimize(PLFvji);
                        PLFiv.minimize(PLFijv); 
                    }   
                }
                shortcuts1[k][i] = PLFvi;
                shortcuts2[k][i] = PLFiv;    
            }
        }    
    }

    void shortcutQueryFunctionUp2Bottom(NodeId s, NodeId d, double t, PLF &costFunction)
    {
        vector<vector<PLF>> sResults1, sResults2, dResults1, dResults2; 
        NodeId lca = qLCA(s,d).first;

        auto itsd = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), d);
        auto itds = find(tnodes[d].ancestors.begin(), tnodes[d].ancestors.end(), s);
        if(itsd != tnodes[s].ancestors.end()) // d is an ancestor of s
        {
            int indexsd = itsd - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[indexsd] == 1)
            {
                //cout << "SHORTCUTS" << endl;
                costFunction = tnodes[s].selectShortcurts1[indexsd];
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                shortcutsForNode(s,sResults1, sResults2);
                costFunction = sResults1[sResults1.size()-1][indexsd];
            }
            
        }
        else if (itds != tnodes[d].ancestors.end())  // s is an ancestor of d
        {
            int indexds = itds - tnodes[d].ancestors.begin();

            if (tnodes[d].shortcurtsPosition2[indexds] == 1)
            {
                //cout << "SHORTCUTS" << endl;
                costFunction = tnodes[d].selectShortcurts2[indexds];
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;
                shortcutsForNode(d,dResults1, dResults2);
                costFunction = dResults2[dResults2.size()-1][indexds];
            }
        }
        else
        {
            int shorctcutFLAG = 1;
            auto it = find(tnodes[s].ancestors.begin(), tnodes[s].ancestors.end(), lca);
            int lcaPosition = it - tnodes[s].ancestors.begin();
            if (tnodes[s].shortcurtsPosition1[lcaPosition] == 0 or tnodes[d].shortcurtsPosition2[lcaPosition] == 0) {shorctcutFLAG =0;}
            for (int j: tnodes[lca].position)
            {
                if(shorctcutFLAG ==0){break;}
                if (tnodes[s].shortcurtsPosition1[j] == 0 or tnodes[d].shortcurtsPosition2[j] == 0)
                {
                    shorctcutFLAG = 0;
                    break;
                }
            }


            if (shorctcutFLAG == 1)
            {
                //cout << "SHORTCUTS" << endl;
                tnodes[d].selectShortcurts2[lcaPosition].compound(tnodes[s].selectShortcurts1[lcaPosition],costFunction,lca);
                for (int j: tnodes[lca].position)
                {
                    PLF PLFsd;
                    tnodes[d].selectShortcurts2[j].compound(tnodes[s].selectShortcurts1[j],PLFsd,tnodes[lca].ancestors[j]);
                    costFunction.minimize(PLFsd);                    
                }
                
            }
            else
            {
                //cout << "WITHOUT SHORTCUTS" << endl;

                shortcutsForNode(s,sResults1, sResults2);
                shortcutsForNode(d,dResults1, dResults2);

                dResults2[dResults2.size()-1][lcaPosition].compound(sResults1[sResults1.size()-1][lcaPosition],costFunction,lca);
                for (int j:tnodes[lca].position)
                {
                    PLF PLFlca;
                    dResults2[dResults2.size()-1][j].compound(sResults1[sResults1.size()-1][j],PLFlca,tnodes[lca].ancestors[j]);
                    costFunction.minimize(PLFlca);
                }          
            }
  
        }
    }

    void shortcutsForNode(NodeId s)
    {
        vector<vector<PLF>> sResults1, sResults2;
        shortcutsForNode(s,sResults1, sResults2);

        vector<NodeId> nodes;
        for (int v = 0; v < tnodes[s].ancestors.size(); v++)
        {
            nodes.push_back(tnodes[s].ancestors[v]);
        }
        nodes.push_back(s);

        for (int k = 0; k<nodes.size();k++)
        {
            NodeId nodev = nodes[k];
            for (int i = 0; i < tnodes[nodev].ancestors.size(); i++)
            {
                tnodes[nodev].shortcurtsPosition1[i] = 1;
                tnodes[nodev].selectShortcurts1[i] = sResults1[k][i];

                tnodes[nodev].shortcurtsPosition2[i] = 1;
                tnodes[nodev].selectShortcurts2[i] = sResults2[k][i];

                //num += 2;
                //cout << nodev+1 << "  " << i << "-th ancestor:  " << tnodes[nodev].selectShortcurts1[i] << "   ||   " << tnodes[nodev].selectShortcurts2[i] << endl;
            }
            
        }  
    }

    void BuildShortcuts(int N)
    {
        int num = 0;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            int height = tnodes[i].height;
            if (height >= height2nodeQ.size()) {
                height2nodeQ.resize(height + 1);
            }
            height2nodeQ[height].push_back(i);            
        }

        while (!height2nodeQ.empty())        
        {
            for (vector<NodeId>::iterator a = height2nodeQ.back().begin(); a!=height2nodeQ.back().end(); a++)
            {
                num ++;
                NodeId vid = (*a);
                auto &node = tnodes[vid];
                cout << "selectShortcutAndBuils: " << vid << " num:" << num<<  endl;
                if (num >= N){break;}
                shortcutsForNode(vid);

            }
            if (num >= N){break;}
            height2nodeQ.pop_back();
        }        
    }
    

    void buildShortcutSetUp2BottomTree(vector<NodeId> &queryNode, int N)
    {
        int shorctnum = 0;
        int verticxnum = 0;
        // initialize all shortcuts 

        cout << "Root: " << root +1 <<endl;
        auto &rootNode = tnodes[root];
        queue<NodeId> Q;
        for (auto nodeid:rootNode.cnodeid)
        {
            tnodes[nodeid].shortcurtsPosition1[0] = 1;
            tnodes[nodeid].selectShortcurts1[0] = tnodes[nodeid].edgesPLF1[0];

            tnodes[nodeid].shortcurtsPosition2[0] = 1;
            tnodes[nodeid].selectShortcurts2[0] = tnodes[nodeid].edgesPLF2[0];

            shorctnum +=2;
            for(auto cnode:tnodes[nodeid].cnodeid){Q.push(cnode);}
        }
        
        while (!Q.empty())
        {
            NodeId v = Q.front();
            Q.pop();
            //cout << "Build shortcuts for:" << v+1 <<"-------------------------------------------------------- " << endl;
            for (int i = 0; i < tnodes[v].ancestors.size(); i++)
            {
                //cout << " -> " << i << "-th ancestor:"<< tnodes[v].ancestors[i]+1  << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<< endl;
                //PLF PLFvi, PLFiv;
                Segment seg(0,0);
                PLF PLFvi, PLFiv;
                PLFvi.f->push_back(seg);
                PLFiv.f->push_back(seg);

                for (int j = 0; j< tnodes[v].edgesId.size(); j++)
                {
                    //cout << tnodes[v].edgesId[j]+1 << "  " << tnodes[v].ancestors[i]+1 << endl;
                    //cout << "through: " << tnodes[v].edgesId[j]+1 << endl;
                    PLF PLFvj, PLFji, PLFvji; // v->j->i
                    PLF PLFij, PLFjv, PLFijv; // i->j->v
                    
                    if(tnodes[v].position[j] == i)
                    {
                        //cout << "case 3: j == i *****************************" << endl; 
                        PLFvji = tnodes[v].edgesPLF1[j];
                        PLFijv = tnodes[v].edgesPLF2[j];

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            // cout << "case 3.1 : i-th shortcut has not been built" << endl;
                            // cout << "shortcut PLFvi: " << PLFvji << endl;
                            // cout << "shortcut PLFiv: " << PLFijv << endl;
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;                            
                        }
                        else
                        {
                            // cout << "case 3.2 : i-th shortcut has been built" << endl;
                            // cout << "shortcut PLFvi before : " << PLFvji << endl;
                            // cout << "shortcut PLFiv before : " << PLFijv << endl;
                            // cout << "PLFvi.minimize(PLFvji);" << endl;
                            PLFvi.minimize(PLFvji);
                            // cout << "PLFiv.minimize(PLFijv);" << endl;
                            PLFiv.minimize(PLFijv);
                            // cout << "new shortcut PLFvi: " << PLFvji << endl;
                            // cout << "new shortcut PLFiv: " << PLFijv << endl;                              
                            // cout << "shortcut PLFvi after : " << PLFvi << endl;
                            // cout << "shortcut PLFiv after : " << PLFiv << endl;                               
                        }
                    }
                    else
                    {
                        if (tnodes[v].position[j] > i)
                        {
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = tnodes[tnodes[v].edgesId[j]].selectShortcurts1[i];

                            PLFij = tnodes[tnodes[v].edgesId[j]].selectShortcurts2[i];
                            PLFjv = tnodes[v].edgesPLF2[j];
                            
                            // cout << "case 1: j>i: v->j->i *****************************" << endl;
                            // cout << "PLFvj : " << PLFvj << endl; 
                            // cout << "PLFji : " << PLFji << endl; 
                            // cout << "PLFij : " << PLFij << endl; 
                            // cout << "PLFjv : " << PLFjv << endl; 
                        }
                        else
                        {
                            //cout << "tnodes[v].position[j] < i" << endl;  
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = tnodes[tnodes[v].ancestors[i]].selectShortcurts2[tnodes[v].position[j]];

                            PLFij = tnodes[tnodes[v].ancestors[i]].selectShortcurts1[tnodes[v].position[j]];
                            PLFjv = tnodes[v].edgesPLF2[j];

                            // cout << "case 2: j < i: v->j->i *****************************" << endl;
                            // cout << "PLFvj : " << PLFvj << endl; 
                            // cout << "PLFji : " << PLFji << endl; 
                            // cout << "PLFij : " << PLFij << endl; 
                            // cout << "PLFjv : " << PLFjv << endl;              
                        }
                        PLFji.compound(PLFvj,PLFvji,tnodes[v].edgesId[j]);
                        PLFjv.compound(PLFij,PLFijv,tnodes[v].edgesId[j]);

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;
                            continue;
                        }
                        PLFvi.minimize(PLFvji);
                        PLFiv.minimize(PLFijv); 
                    }   
                    

                }
                // cout << PLFvi << endl;
                // cout << PLFiv << endl;
                tnodes[v].shortcurtsPosition1[i] = 1;
                tnodes[v].selectShortcurts1[i] = PLFvi;
                tnodes[v].shortcurtsPosition2[i] = 1;
                tnodes[v].selectShortcurts2[i] = PLFiv;                
                shorctnum+=2; 
            }
            verticxnum ++;
            for(auto cnode:tnodes[v].cnodeid){Q.push(cnode);}
            cout << "build shortcut for node: " << v << " num:" << shorctnum << " verticxnum "  << verticxnum<< endl;
            // if (verticxnum >= N)
            // {
            //     break;
            // }
            queryNode.push_back(v);
            
        }
        // cout << "build shortcut: " << shorctnum << endl;
    }

    void WriteShortcuts()
    {
        ofstream record_shortcut(record_shortcut_path);
        cout << "Root: " << root +1 <<endl;
        auto &rootNode = tnodes[root];
        queue<NodeId> Q;
        for (auto nodeid:rootNode.cnodeid)
        {
            tnodes[nodeid].shortcurtsPosition1[0] = 1;
            tnodes[nodeid].selectShortcurts1[0] = tnodes[nodeid].edgesPLF1[0];

            tnodes[nodeid].shortcurtsPosition2[0] = 1;
            tnodes[nodeid].selectShortcurts2[0] = tnodes[nodeid].edgesPLF2[0];

            for(auto cnode:tnodes[nodeid].cnodeid){Q.push(cnode);}
        }
        
        while (!Q.empty())
        {
            NodeId v = Q.front();
            Q.pop();
            //cout << "Build shortcuts for:" << v+1 <<"-------------------------------------------------------- " << endl;
            record_shortcut << v << " " <<tnodes[v].height << " " << tnodes[v].ancestors.size() << "\n";
            for (int i = 0; i< tnodes[v].ancestors.size(); i++)
            {
                record_shortcut << tnodes[v].ancestors[i] << " " << tnodes[v].selectShortcurts1[i].f->size() << "\n";
                auto seg = tnodes[v].selectShortcurts1[i].f->begin();
                while (seg!= tnodes[v].selectShortcurts1[i].f->end())
                {
                    record_shortcut << seg->t << " " << seg->w <<"\n";
                    seg++;
                }

            }
            
            for(auto cnode:tnodes[v].cnodeid){Q.push(cnode);}
            cout << "Record shortcut for node: " << v <<  endl;
            // if (verticxnum >= N)
            // {
            //     break;
            // }
            
        }
        record_shortcut.close();

    }

    void BuildandWriteShorcuts()
    {
        int shorctnum = 0;
        int verticxnum = 0;
        // initialize all shortcuts 
        
        ofstream record_shortcut(record_shortcut_path);
        cout << "Root: " << root +1 <<endl;
        auto &rootNode = tnodes[root];
        queue<NodeId> Q;
        for (auto nodeid:rootNode.cnodeid)
        {
            tnodes[nodeid].shortcurtsPosition1[0] = 1;
            tnodes[nodeid].selectShortcurts1[0] = tnodes[nodeid].edgesPLF1[0];

            tnodes[nodeid].shortcurtsPosition2[0] = 1;
            tnodes[nodeid].selectShortcurts2[0] = tnodes[nodeid].edgesPLF2[0];

            shorctnum +=2;
            for(auto cnode:tnodes[nodeid].cnodeid){Q.push(cnode);}
        }
        
        while (!Q.empty())
        {
            NodeId v = Q.front();
            Q.pop();
            
            // if (tnodes[v].height >=10)
            // {
            //     break;
            // }
            

            record_shortcut << v << " " <<tnodes[v].height << " " << tnodes[v].ancestors.size() << "\n";
            for (int i = 0; i < tnodes[v].ancestors.size(); i++)
            {
                Segment seg(0,0);
                PLF PLFvi, PLFiv;
                PLFvi.f->push_back(seg);
                PLFiv.f->push_back(seg);

                for (int j = 0; j< tnodes[v].edgesId.size(); j++)
                {
                    PLF PLFvj, PLFji, PLFvji; // v->j->i
                    PLF PLFij, PLFjv, PLFijv; // i->j->v
                    
                    if(tnodes[v].position[j] == i)
                    {
                        //cout << "case 3: j == i *****************************" << endl; 
                        PLFvji = tnodes[v].edgesPLF1[j];
                        PLFijv = tnodes[v].edgesPLF2[j];

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;                            
                        }
                        else
                        {
                            PLFvi.minimize(PLFvji);
                            // cout << "PLFiv.minimize(PLFijv);" << endl;
                            PLFiv.minimize(PLFijv);                            
                        }
                    }
                    else
                    {
                        if (tnodes[v].position[j] > i)
                        {
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = tnodes[tnodes[v].edgesId[j]].selectShortcurts1[i];

                            PLFij = tnodes[tnodes[v].edgesId[j]].selectShortcurts2[i];
                            PLFjv = tnodes[v].edgesPLF2[j];
                            
                        }
                        else
                        {
                            PLFvj = tnodes[v].edgesPLF1[j];
                            PLFji = tnodes[tnodes[v].ancestors[i]].selectShortcurts2[tnodes[v].position[j]];

                            PLFij = tnodes[tnodes[v].ancestors[i]].selectShortcurts1[tnodes[v].position[j]];
                            PLFjv = tnodes[v].edgesPLF2[j];            
                        }
                        PLFji.compound(PLFvj,PLFvji,tnodes[v].edgesId[j]);
                        PLFjv.compound(PLFij,PLFijv,tnodes[v].edgesId[j]);

                        if (PLFvi.dpt2arr(0) == 0)
                        {
                            PLFvi = PLFvji;
                            PLFiv = PLFijv;
                            continue;
                        }
                        PLFvi.minimize(PLFvji);
                        PLFiv.minimize(PLFijv); 
                    }   
                    

                }
                tnodes[v].shortcurtsPosition1[i] = 1;
                tnodes[v].selectShortcurts1[i] = PLFvi;
                tnodes[v].shortcurtsPosition2[i] = 1;
                tnodes[v].selectShortcurts2[i] = PLFiv;                
                shorctnum+=2; 
            }
            
            for (int i = 0; i< tnodes[v].ancestors.size(); i++)
            {
                record_shortcut << tnodes[v].ancestors[i] << " " << tnodes[v].selectShortcurts1[i].f->size() << "\n";
                //cout << "shortcut size: " << tnodes[v].selectShortcurts1[i].f->size() << endl;
                auto seg = tnodes[v].selectShortcurts1[i].f->begin();
                while (seg!= tnodes[v].selectShortcurts1[i].f->end())
                {
                    record_shortcut << seg->t << " " << seg->w <<"\n";
                    seg++;
                }

            }            
            
            verticxnum ++;
            for(auto cnode:tnodes[v].cnodeid){Q.push(cnode);}
            cout << "build shortcut for node: " << v << "tree height: " << tnodes[v].height <<  " num:" << shorctnum << " verticxnum "  << verticxnum<< endl;
            
        }

        record_shortcut.close();        
    }

    void Init() 
    {

        load_graph_zy();

        vector<int> dgr;

        // for (int i = 0; i < graph.size(); i++)
        // {
        //     dgr.push_back(graph[i].size());
        // }
        // auto maxPosition = max_element(dgr.begin(), dgr.end());
        // auto minPosition = min_element(dgr.begin(), dgr.end());
        // cout << *maxPosition << " at the position of " << maxPosition - dgr.begin() << endl;
        // cout << *minPosition << " at the position of " << minPosition - dgr.begin() << endl;

        //ofstream write_shortcut("result.txt");

        clock_t tbegin = clock();

        creat_index();
        
        //report();
        //restore_index();
        //ShowTreeStructure();

        //return;


        clock_t tend = clock();
        cout << "Finish indexing time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;
        //ofstream index_result(result_path);
        //index_result << "Finish indexing time cost: " <<(tend - tbegin) / CLOCKS_PER_SEC << "\n"; 

        cout << endl;
        cout << endl;


        int h = 0,w = 0;
        for (int i = 0; i < tnodes.size(); i++)
        {
            if (tnodes[i].height > h)
            {
                h = tnodes[i].height;
            }
            if (tnodes[i].edgesId.size() > w)
            {
                w = tnodes[i].edgesId.size();
            }            
        }

        cout << "tree height: " << h << "  tree width: " << w << endl;

        ///////////////////////////// Parameters setting
        int N = 25000000000000000000000000;

    
        clock_t tbegin3 = clock();
        //buildShortcutsSet();
        //buildShortcutSetUp2Bottom();
        vector<NodeId> queryNode;
        //buildShortcutSetUp2BottomTree(queryNode, N);
        BuildandWriteShorcuts();
        cout << "Finish building all shortcuts" << endl;
        //WriteShortcuts();
        // for (unsigned int i = 0; i < tnodes.size(); i++)
        // {
        //     int leafnum = 0;
        //     if (tnodes[i].cnodeid.size() == 0)
        //     {
        //         leafnum ++;
        //         cout << "1: ";
        //         cout << tnodes[i].selectShortcurts1[tnodes[i].selectShortcurts1.size()-1] << endl;
        //         cout << "2: " ;
        //         cout << tnodes[i].selectShortcurts1[(tnodes[i].selectShortcurts1.size()-1)/2] << endl;
        //         cout << "3: " ;
        //         cout << tnodes[i].selectShortcurts1[0] << endl;
        //         cout << " -------------------------------------------------- " << endl;
        //         cout << endl;
        //     }
        //     if (leafnum >= 50)
        //     {
        //         break;
        //     }            
        // }
    }

};




#endif //TREE_DECOMP_H



