#include <stdio.h>

#include <Rembedded.h>
#include <Rcpp.h>

using namespace Rcpp;

// have to replicate this here, as RcppExports isn’t a header
List knn_asym(const NumericMatrix data, const size_t k, const std::string distance);

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

void main_impl() {
	setup_packages();
	NumericMatrix guo = load_guo();
	
	List knns = knn_asym(guo, 5, "euclidean");
	NumericMatrix dists = knns["dist"];
	std::cout << dists << std::endl;
}

int main(int argc, char* argv[]) {
	char* argv_r[] = {argv[0], (char*)"--vanilla", (char*)"--slave"};
	Rf_initEmbeddedR(sizeof(argv_r) / sizeof(argv_r[0]), argv_r);
	main_impl();
	Rf_endEmbeddedR(0);
	return 0;
}
