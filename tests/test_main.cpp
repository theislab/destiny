#include <iostream>
#include <fstream>

#include <Rembedded.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/SparseExtra>

#include "../src/exports.h"

using namespace Rcpp;


void setup_packages() {
	Environment base("package:base");
	Function library = base["library"];
	
	library("Rcpp");
	library("Matrix");
}

NumericMatrix exprs(RObject exprset) {
	Environment assay_data = exprset.attr("assayData");
	return assay_data["exprs"];
}

NumericMatrix load_guo() {
	Environment global = Environment::global_env();
	
	Environment utils("package:utils");
	Function load_data = utils["data"];
	
	load_data("guo_norm", Named("package", "destiny"));
	return exprs(global["guo_norm"]);
}

void save_tsv(NumericMatrix m, std::string filename) {
	Eigen::IOFormat TSVFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n", "", "", "", "");
	Eigen::MatrixXd me(Rcpp::as<Eigen::MatrixXd>(m));
	
	std::ofstream f;
	f.open(filename);
	f << me.format(TSVFmt) << std::endl;
	f.close();
}

void main_impl() {
	setup_packages();
	NumericMatrix guo = load_guo();
	
	List knns = knn_asym(guo, 5, "euclidean");
	NumericMatrix index = knns["index"];
	NumericMatrix dist = knns["dist"];
	Eigen::SparseMatrix<double> dist_mat = knns["dist_mat"];
	
	save_tsv(index, "index.tsv");
	save_tsv(dist, "dist.tsv");
	Eigen::saveMarket(dist_mat, "dist_mat.mm");
}

int main(int argc, char* argv[]) {
	char* argv_r[] = {argv[0], (char*)"--vanilla", (char*)"--slave"};
	Rf_initEmbeddedR(sizeof(argv_r) / sizeof(argv_r[0]), argv_r);
	main_impl();
	Rf_endEmbeddedR(0);
	return 0;
}
