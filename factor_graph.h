#ifndef FACTOR_GRAPH_H
#define FACTOR_GRAPH_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct Conn {
    std::vector<Conn *> list;
    int vertex;
    float node_val;
};

class Factor_graph {

public:
    Factor_graph();
    ~Factor_graph();
    virtual void create_adj_mat() = 0;
    virtual void create_list_from_mat() = 0;
    void mat_from_file(std::string, int, int);
    void setEdge(int, int, int);
    int getEdge(int, int);
    void setMatSize(int, int);
    void setRowsandCols(int, int);
    int getNumRows();
    int getNumCols();
    std::vector<std::vector<int> >& getAdjMat();
    void setAdjMat(std::vector<std::vector<int> > &);
    
private:
    std::vector<std::vector<int> > adj_mat;
    int rows, cols;
};

#endif