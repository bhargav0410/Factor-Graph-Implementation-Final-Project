#include "factor_graph.h"

Factor_graph::Factor_graph() {}

Factor_graph::~Factor_graph() {}

//gets adjacency matrix from file
void Factor_graph::mat_from_file(std::string file, int numRows, int numCols) {
    setMatSize(numRows, numCols);
    std::ifstream infile;
    int mat_size;
    infile.open(file.c_str(), std::ifstream::binary);
    if (infile.is_open()) {
        infile.seekg(0, infile.end);
        mat_size = infile.tellg()/sizeof(int);
        infile.seekg(0, infile.beg);
        if (mat_size > numCols*numRows) {
            std::cerr << "Re-check matrix size...\n";
        } else {
            for (int i = 0; i < numRows; i++) {
                infile.read((char *)&adj_mat[i][0], numCols*sizeof(int));
            }
        }
    } else {
        std::cerr << "File does not exist...\n";
    }
}

void Factor_graph::setEdge(int row, int col, int val) {
    if (row <= rows and col <= cols) {
        adj_mat[row][col] = val;
    } else {
        std::cerr << "Out of bounds...\n";
    }
}

int Factor_graph::getEdge(int row, int col) {
    if (row <= rows and col <= cols) {
        return adj_mat[row][col];
    } else {
        std::cerr << "Out of bounds...\n";
        return 0;
    }
}

void Factor_graph::setMatSize(int numRows, int numCols) {
    setRowsandCols(numRows, numCols);
    adj_mat.resize(rows);
    for (int i = 0; i < rows; i++) {
        adj_mat[i].resize(cols);
    }
}

void Factor_graph::setRowsandCols(int numRows = 0, int numCols = 0) {
    if (numRows > 0) {
        rows = numRows;
    }
    if (numCols > 0) {
        cols = numCols;
    }
}

int Factor_graph::getNumRows() {
    return rows;
}

int Factor_graph::getNumCols() {
    return cols;
}

std::vector<std::vector<int> >& Factor_graph::getAdjMat() {
    return adj_mat;
}

void Factor_graph::setAdjMat(std::vector<std::vector<int> > &vec) {
    adj_mat = vec;
}
