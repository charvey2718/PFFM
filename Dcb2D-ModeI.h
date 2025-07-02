#ifndef DCB2D_H
#define DCB2D_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <stdio.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

std::vector<double> linspace(double start, double end, double initial_size, double growth_rate = 1, double tol = 1e-9);
std::vector<double> linspace2(double start, double end, double initial_size, double growth_rate, double max_size, double tol = 1e-9);
void insert_midpoints(std::vector<double>& vec);
void meshgrid(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y);
void create_nodes_matrix(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& nodes);
void create_elements_matrix(const std::vector<std::vector<double>>& X, const std::vector<std::vector<double>>& Y, bool isQuadratic, std::vector<std::vector<int>>& elements);
void add_crack(std::vector<std::vector<double>>& nodes, std::vector<std::vector<int>>& elements, double xct, double yc, double tol = 1e-12);

void Get2DQuadrature(int degree, Eigen::VectorXd& xis, Eigen::VectorXd& etas, Eigen::VectorXd& ws);

Eigen::MatrixXd LinearSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double E, double nu, double t,Eigen::VectorXd u, Eigen::VectorXd s, bool plane_stress, int quadratureDegree);
Eigen::MatrixXd QuadraticSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double E, double nu, double t, Eigen::VectorXd u, Eigen::VectorXd s, bool plane_stress, int quadratureDegree);
void LinearPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Gc, double ls, Eigen::VectorXd s, std::vector<double> H, int quadratureDegree, Eigen::MatrixXd& k, Eigen::VectorXd& f);
void QuadraticPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double Gc, double ls, Eigen::VectorXd s, std::vector<double> H, int quadratureDegree, Eigen::MatrixXd& k, Eigen::VectorXd& f);

void AssembleLinearSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m);
void AssembleLinearPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m);
void AssembleLinearPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m);

void AssembleQuadraticSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r);
void AssembleQuadraticPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r);
void AssembleQuadraticPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m, int n, int o, int p, int q, int r);

void AssembleAll(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, const Eigen::VectorXd& s, const std::vector<std::vector<double>>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double Gc_eff, double Gc_bulk, double t, double h2, double hi, double ls, Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, bool isPlaneStress, bool isQuadratic, int quadratureDegree);

void UpdateH(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, std::vector<std::vector<double>>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double h2, double hi, double ls, bool isPlaneStress, bool isQuadratic, int quadratureDegree);
std::vector<double> GetLinearElementTensileElasticStrainEnergy(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double E, double nu, Eigen::VectorXd u, bool plane_stress, int quadratureDegree);
std::vector<double> GetQuadraticElementTensileElasticStrainEnergy(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double E, double nu, Eigen::VectorXd u, bool plane_stress, int quadratureDegree);
void SaveHField(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const std::vector<std::vector<double>>& H, bool isQuadratic, int quadratureDegree, std::string filename, int prec);

double GetGcFromComplianceCalibration(const std::vector<double>& D, const std::vector<double>& P, const std::vector<double>& a, double a_min);
double GetGcFromdPidA(const std::vector<double>& D, const std::vector<double>& P, const std::vector<double>& a, const std::vector<double>& U, double a_min);

#endif
