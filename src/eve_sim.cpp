//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#include <limits>
#include <algorithm>
#include <numeric>
#include <string>
#include "evesim/rutils.h"
#include "evesim/phylo.h"
#include "evesim/tree_metric.h"
#include "evesim/sim_table.h"


using namespace rutils;
using namespace tres_sim;


namespace {

  constexpr const char* Xtree_tag = "tres_sim::Xtree_tag";
  constexpr const char* SimTable_tag = "tres_sim::SimTable_tag";


  // descending LR
  ltable_t L2_ltable(const Rcpp::NumericMatrix& LR) {
    if ((4 < LR.ncol()) || (2 > LR.nrow())) {
      Rcpp::stop("Matrix \"LR\" is not an conformant Ltable");
    }
    return ltable_view<true>(LR, LR(0, 0)).to_ltable();
  }


  // autoselect descending/ascending LR
  ltable_t L2_ltable(const Rcpp::NumericMatrix& LR, Rcpp::Nullable<double> age) {
    if ((4 < LR.ncol()) || (2 > LR.nrow())) {
      Rcpp::stop("Matrix \"LR\" is not an conformant Ltable");
    }
    if (0.0 == LR(0,0)) { // ascending
      if (age.isNull()) {
        Rcpp::stop("Argument \"age\" required for ascending Ltable");
      }
      return ltable_view<false>(LR, Rcpp::as<double>(age)).to_ltable();
    }
    // descending
    const auto eff_age = std::min<double>(LR(0,0), age.isNotNull() ? Rcpp::as<double>(age) : std::numeric_limits<int>::max());
    return ltable_view<true>(LR, eff_age).to_ltable();
  }

}


// [[Rcpp::export("Ltable.phylo", rng = false)]]
Rcpp::List Ltable_phylo(const Rcpp::NumericMatrix& LR,
                        bool drop_extinct,
                        Rcpp::Nullable<double> age = R_NilValue) {
  auto ltable = L2_ltable(LR, age);
  return phylo_t(tree_t::from(ltable, drop_extinct)).unwrap();
}


// [[Rcpp::export("Ltable.prune", rng = false)]]
Rcpp::NumericMatrix Ltable_prune(Rcpp::NumericMatrix LR, Rcpp::Nullable<double> age = R_NilValue) {
  if ((4 < LR.ncol()) || (2 > LR.nrow())) {
    Rcpp::stop("Matrix \"LR\" is not an conformant Ltable");
  }
  if (0.0 == LR(0,0)) { // ascending
    if (age.isNull()) {
      Rcpp::stop("Argument \"age\" required for ascending Ltable");
    }
    return ltable_view<false>(LR, Rcpp::as<double>(age)).to_matrix();
  }
  // descending
  const auto eff_age = std::min<double>(LR(0,0), age.isNotNull() ? Rcpp::as<double>(age) : std::numeric_limits<int>::max());
  return ltable_view<true>(LR, eff_age).to_matrix();
}


// [[Rcpp::export("Ltable.legacy.ascending", rng = false)]]
Rcpp::NumericMatrix Ltable_legacy_ascending(Rcpp::NumericMatrix LR, Rcpp::Nullable<double> age = R_NilValue) {
  if (0.0 != LR(0,0)) {
    // descending source
    double t = LR(0,0);
    auto L = ltable_view<true>(LR, t).to_matrix();
    for (int i = 0; i < L.rows(); ++i) {
      L(i, 0) = t - LR(i, 0);
      L(i, 3) = (L(i, 3) == -1.0) ? -1.0 : t - L(i, 3);
    }
    L(1,1) = -1;  // reintroduce legacy quirk
    return L;
  }
  // ascending source
  auto L = ltable_view<false>(LR, Rcpp::as<double>(age)).to_matrix();
  L(1,1) = -1;  // reintroduce legacy quirk
  return L;
}


// [[Rcpp::export("Ltable.legacy.descending", rng = false)]]
Rcpp::NumericMatrix Ltable_legacy_descending(Rcpp::NumericMatrix LR, Rcpp::Nullable<double> age = R_NilValue) {
  if (0.0 != LR(0,0)) {
    // descending source
    auto L = ltable_view<true>(LR, LR(0,0)).to_matrix();
    L(1,1) = -1;  // reintroduce legacy quirk
    return L;
  }
  // ascending source
  const auto t = Rcpp::as<double>(age);
  auto L = ltable_view<false>(LR, t).to_matrix();
  for (int i = 0; i < L.rows(); ++i) {
    L(i, 0) = t - LR(i, 0);
    L(i, 3) = (L(i, 3) == -1.0) ? -1.0 : t - L(i, 3);
  }
  L(1,1) = -1;  // reintroduce legacy quirk
  return L;
}


// [[Rcpp::export("Ltable.cophenetic", rng = false)]]
Rcpp::NumericMatrix Ltable_cophenetic(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto tree = tree_t::from(L2_ltable(LR, age), drop_extinct);
  return tree_metric::metric(tree_metric::cophenetic{}, tree);
}


// [[Rcpp::export("Ltable.ed", rng = false)]]
Rcpp::NumericVector Ltable_ed(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto tree = tree_t::from(L2_ltable(LR, age), drop_extinct);
  return tree_metric::metric(tree_metric::ed{}, tree);
}


// [[Rcpp::export("Ltable.nnd", rng = false)]]
Rcpp::NumericVector Ltable_nnd(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto tree = tree_t::from(L2_ltable(LR, age), drop_extinct);
  return tree_metric::metric(tree_metric::nnd{}, tree);
}


// [[Rcpp::export("Ltable.pd", rng = false)]]
double Ltable_pd(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto tree = tree_t::from(L2_ltable(LR, age), drop_extinct);
  return tree_metric::metric(tree_metric::pd{}, tree);
}


// [[Rcpp::export("Ltable.mpd", rng = false)]]
double Ltable_mpd(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto tree = tree_t::from(L2_ltable(LR, age), drop_extinct);
  return tree_metric::metric(tree_metric::mpd{}, tree);
}


// [[Rcpp::export("Ltable.tree", rng = false)]]
SEXP Ltable_tree(Rcpp::NumericMatrix LR, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  auto ltable = L2_ltable(LR, age);
  auto ptree = new tree_t{tree_t::from(ltable, drop_extinct)};
  return tagged_xptr<tree_t>(ptree, Xtree_tag);
}


// [[Rcpp::export("Xtree", rng = false)]]
SEXP Xtree(Rcpp::RObject Robj, bool drop_extinct, Rcpp::Nullable<double> age = R_NilValue) {
  if (Rcpp::is<Rcpp::NumericMatrix>(Robj)) {
    // Ltable
    return Ltable_tree(Rcpp::as<Rcpp::NumericMatrix>(Robj), drop_extinct, age);
  }
  if (!age.isNull()) {
    Rcpp::warning("Xtree(): argument \"age\" ignored");
  }
  if (Robj.inherits("phylo")) {
    // ape::phylo
    auto phy = Rcpp::clone(Rcpp::as<Rcpp::List>(Robj));
    if (0 == Rcpp::as<int>(phy["Nnode"])) Rcpp::stop("Xtree(): insufficient nodes in phylo object");
    if (std::string(phy.attr("order")) != "cladewise") {
      Rcpp::Function reorder("reorder", Rcpp::Environment::namespace_env("ape"));
      phy = reorder(phy);
    }
    if (!Rcpp::is<Rcpp::IntegerVector>(phy["tip.label"])) {
      auto nnode = Rcpp::as<int>(phy["Nnode"]);
      auto labels = Rcpp::IntegerVector(nnode + 1);
      std::iota(labels.begin(), labels.end(), 0);
      phy["tip.label"] = labels;
    }
    auto ptree = new tree_t(tree_t::from(phylo_t(std::move(phy))));
    return tagged_xptr<tree_t>(ptree, Xtree_tag);
  }
  // Xtree, clone
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  auto pclone = new tree_t(*ptr);
  return tagged_xptr<tree_t>(pclone, Xtree_tag);
}


// [[Rcpp::export("Xtree.is", rng = false)]]
bool Xtree_is(Rcpp::RObject Robj) {
  return is_tagged_xptr<tree_t>(Robj, Xtree_tag);
}


// [[Rcpp::export("Xtree.tips", rng = false)]]
int Xtree_tips(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return ptr->tips();
}


// [[Rcpp::export("Xtree.nnode", rng = false)]]
int Xtree_nnode(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return ptr->nnode();
}


// [[Rcpp::export("Xtree.drop.extinct", rng = false)]]
SEXP Xtree_drop_extinct(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  auto tree = drop_tips_if(*ptr, [](const auto& node) { return node.t != 0.0; });
  auto ptree = new tree_t(std::move(tree));
  return tagged_xptr<tree_t>(ptree, Xtree_tag);
}


// [[Rcpp::export("Xtree.tip.label", rng = false)]]
Rcpp::IntegerVector Xtree_tip_label(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  auto tip_label = Rcpp::IntegerVector(ptr->tips());
  auto& nodes = ptr->nodes;
  for (int i = 0; i < tip_label.size(); ++i) {
    tip_label[i] = std::abs(nodes[i].label);
  }
  return tip_label;
}


// [[Rcpp::export("Xtree.phylo", rng = false)]]
Rcpp::List Xtree_phylo(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return phylo_t(*ptr).unwrap();
}


// [[Rcpp::export("Xtree.cophenetic", rng = false)]]
Rcpp::NumericMatrix Xtree_cophenetic(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return tree_metric::metric(tree_metric::cophenetic{}, *ptr);
}


// [[Rcpp::export("Xtree.ed", rng = false)]]
Rcpp::NumericVector Xtree_ed(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return tree_metric::metric(tree_metric::ed{}, *ptr);
}


// [[Rcpp::export("Xtree.nnd", rng = false)]]
Rcpp::NumericVector Xtree_nnd(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return tree_metric::metric(tree_metric::nnd{}, *ptr);
}


// [[Rcpp::export("Xtree.pd", rng = false)]]
double Xtree_pd(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return tree_metric::metric(tree_metric::pd{}, *ptr);
}


// [[Rcpp::export("Xtree.mpd", rng = false)]]
double Xtree_mpd(SEXP Robj) {
  auto* ptr = tagged_xptr<tree_t>(Robj, Xtree_tag).get();
  return tree_metric::metric(tree_metric::mpd{}, *ptr);
}


// [[Rcpp::export("SimTable", rng = false)]]
SEXP SimTable(Rcpp::Nullable<Rcpp::RObject> Rrhs = R_NilValue) {
  auto Robj = Rcpp::as<Rcpp::RObject>(Rrhs);
  if (Rcpp::is<Rcpp::NumericMatrix>(Robj)) {
    auto L = Rcpp::as<Rcpp::NumericMatrix>(Robj);
    auto age = L(0,0);
    if (age == 0.0) {
      Rcpp::stop("SimTable(): descending Ltable expected");
    }
    auto ltable = L2_ltable(L);
    return tagged_xptr<sim_table_t>(new sim_table_t(std::move(ltable)), SimTable_tag);
  }
  if (is_tagged_xptr<sim_table_t>(Robj, SimTable_tag)) {
    auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
    return tagged_xptr<sim_table_t>(new sim_table_t(*ptr), SimTable_tag);
  }
  if (Rcpp::is<double>(Rrhs)) {
    return tagged_xptr<sim_table_t>(new sim_table_t(Rcpp::as<double>(Rrhs)), SimTable_tag);
  }
  Rcpp::stop("invalid argument to \"SimTable\"");
}


// [[Rcpp::export("SimTable.is", rng = false)]]
bool SimTable_is(Rcpp::RObject Robj) {
  return is_tagged_xptr<sim_table_t>(Robj, SimTable_tag);
}


// [[Rcpp::export("SimTable.age", rng = false)]]
double SimTable_age(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  return ptr->age();
}


// [[Rcpp::export("SimTable.nspecie", rng = false)]]
int SimTable_nspecie(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  return ptr->nspecie();
}


// [[Rcpp::export("SimTable.nclade_specie", rng = false)]]
Rcpp::IntegerVector SimTable_nclade_specie(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto cs = ptr->nclade_specie();
  return Rcpp::IntegerVector(cs.begin(), cs.end());
}


// [[Rcpp::export("SimTable.size", rng = false)]]
int SimTable_size(Rcpp::RObject Robj) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  return ptr->size();
}


//' Convert SimTable object to Legacy L-table Format
//'
//' This function extracts the simulation results from a SimTable object and converts them into
//' a legacy L-table format, represented as an R numeric matrix.
//'
//' The L-table is a common representation of evolutionary simulations, where each row corresponds to a lineage
//' with its birth, ancestor, label, and death times. This function also preserves a known flaw in the legacy L-table
//' format where the second column's first element is manually set to -1.
//'
//' @param Robj A tagged external pointer to a `sim_table_t` object, which contains the SimTable simulation results.
//'
//' @return A \code{NumericMatrix} with four columns:
//' \describe{
//'   \item{Column 1}{The birth time (t) of the lineage.}
//'   \item{Column 2}{The ancestor of the lineage.}
//'   \item{Column 3}{The label assigned to the lineage.}
//'   \item{Column 4}{The death time of the lineage, or -1 if the lineage is still extant.}
//' }
//'
//' @note For backward compatibility, the second element of the first row in the matrix (Column 2, Row 1) is set to -1 to reintroduce a known flaw in the legacy L-table format used in the DDD package.
//'
//' @examples
//' # Example usage
//' pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
//' sim_table <- evesim::edd_sim(pars = pars, age = 20, metric = "nnd", offset = "none")
//' ltable <- SimTable.ltable(sim_table$sim)
//'
//' @export
//' @useDynLib evesim
// [[Rcpp::export("SimTable.ltable", rng = false)]]
Rcpp::NumericMatrix SimTable_ltable(SEXP Robj) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  const auto& ltable = ptr->ltable();
  auto LR = Rcpp::NumericMatrix(ltable.size(), 4);
  for (int i = 0; i < ltable.size(); ++i) {
    auto LRi = LR.row(i);
    LRi[0] = ltable[i].t;
    LRi[1] = (ltable[i].label < 0) ? -(ltable[i].ancestor + 1): (ltable[i].ancestor + 1);
    LRi[2] = ltable[i].label;
    LRi[3] = (ltable[i].death == 0.0) ? -1.0 : ltable[i].death;
  }
  // reintroduce legacy Ltable flaw:
  //LR(1,1) = -1;
  return LR;
}


//' Convert SimTable object to Phylogenetic Tree
//'
//' This function extracts a phylogenetic tree from a given SimTable object and returns it in a
//' phylo-compatible format, suitable for downstream analysis in R. It can optionally remove extinct lineages
//' from the tree.
//'
//' @param Robj A tagged external pointer to a `sim_table_t` object, which contains the SimTable simulation results.
//' @param drop_extinct A logical flag indicating whether extinct lineages should be removed from the resulting phylogenetic tree.
//' \code{TRUE} will drop extinct lineages, and \code{FALSE} will return the full tree including extinct lineages.
//'
//' @return A phylo object (from the \code{ape} package) representing the phylogenetic tree extracted from the simulation table.
//'
//' @examples
//' # Example usage
//' pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
//' sim_table <- evesim::edd_sim(pars = pars, age = 20, metric = "nnd", offset = "none")
//' phylo_tree <- SimTable.phylo(sim_table$sim, drop_extinct = TRUE)
//'
//' @seealso \code{\link[ape]{phylo}} for more information about phylogenetic trees in R.
//'
//' @export
//' @useDynLib evesim
// [[Rcpp::export("SimTable.phylo", rng = false)]]
Rcpp::List SimTable_phylo(SEXP Robj, bool drop_extinct) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  return drop_extinct
    ? phylo_t(ptr->tree()).unwrap()
    : phylo_t(ptr->full_tree()).unwrap();
}


// [[Rcpp::export("SimTable.tree", rng = false)]]
SEXP SimTable_tree(SEXP Robj, bool drop_extinct) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  return drop_extinct
    ? tagged_xptr<tree_t>(new tree_t(ptr->tree()), Xtree_tag)
    : tagged_xptr<tree_t>(new tree_t(ptr->full_tree()), Xtree_tag);
}


// [[Rcpp::export("SimTable.speciation", rng = false)]]
Rcpp::IntegerVector SimTable_speciation(SEXP Robj, int specie, double t) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  --specie;
  if ((specie < 0 ) || (specie > ptr->nspecie())) {
    Rcpp::stop("SimTable.speciation: \"specie\" not in ltable");
  }
  ptr->speciation(specie, t);
  auto ncs = ptr->nclade_specie();
  return Rcpp::IntegerVector(ncs.begin(), ncs.end());
}


// [[Rcpp::export("SimTable.extinction", rng = false)]]
Rcpp::IntegerVector SimTable_extinction(SEXP Robj, int specie, double t) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  --specie;
  if ((specie < 0 ) || (specie > ptr->nspecie())) {
    Rcpp::stop("SimTable.extinction: \"specie\" not in ltable");
  }
  ptr->extinction(specie, t);
  auto ncs = ptr->nclade_specie();
  return Rcpp::IntegerVector(ncs.begin(), ncs.end());
}


// [[Rcpp::export("SimTable.tip.label", rng = false)]]
Rcpp::IntegerVector SimTable_tip_label(SEXP Robj, bool drop_extinct) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  if (!drop_extinct) {
    auto tree = ptr->full_tree();
    auto ret = Rcpp::IntegerVector(tree.tips());
    tree.tip_label(ret.begin());
    return ret;
  }
  auto ret = Rcpp::IntegerVector(ptr->nspecie());
  ptr->specie_label(ret.begin());
  return ret;
}


// [[Rcpp::export("SimTable.cophenetic", rng = false)]]
Rcpp::NumericMatrix SimTable_cophenetic(SEXP Robj, Rcpp::Nullable<double> t = R_NilValue) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto res = tree_metric::metric(tree_metric::cophenetic{}, ptr->tree());
  if (t.isNotNull()) {
    const double ofs = 2.0 * (ptr->age() - Rcpp::as<double>(t));
    const int tips = ptr->nspecie();
    for (int i = 0; i < tips; ++i) {
      for (int j = 0; j < tips; ++j) {
        if (i != j) res(i, j) -= ofs;
      }
    }
  }
  return res;
}


// [[Rcpp::export("SimTable.ed", rng = false)]]
Rcpp::NumericVector SimTable_ed(SEXP Robj, Rcpp::Nullable<double> t = R_NilValue) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto res =tree_metric::metric(tree_metric::ed{}, ptr->tree());
  if (t.isNotNull()) {
    const double ofs = 2.0 * (ptr->age() - Rcpp::as<double>(t));
    for (auto& r : res) { r -= ofs; }
  }
  return res;
}


// [[Rcpp::export("SimTable.nnd", rng = false)]]
Rcpp::NumericVector SimTable_nnd(SEXP Robj, Rcpp::Nullable<double> t = R_NilValue) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto res = tree_metric::metric(tree_metric::nnd{}, ptr->tree());
  if (t.isNotNull()) {
    const double ofs = 2.0 * (ptr->age() - Rcpp::as<double>(t));
    for (auto& r : res) { r -= ofs; }
  }
  return res;
}


// [[Rcpp::export("SimTable.pd", rng = false)]]
double SimTable_pd(SEXP Robj, Rcpp::Nullable<double> t = R_NilValue) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto res = tree_metric::metric(tree_metric::pd{}, ptr->tree());
  if (t.isNotNull()) {
    res -= ptr->nspecie() * (ptr->age() - Rcpp::as<double>(t));
  }
  return res;
}


// [[Rcpp::export("SimTable.mpd", rng = false)]]
double SimTable_mpd(SEXP Robj, Rcpp::Nullable<double> t = R_NilValue) {
  auto* ptr = tagged_xptr<sim_table_t>(Robj, SimTable_tag).get();
  auto res =  tree_metric::metric(tree_metric::mpd{}, ptr->tree());
  if (t.isNotNull()) {
    const int tips = ptr->nspecie();
    res -= 2.0 * (ptr->age() - Rcpp::as<double>(t)) / (tips * (tips - 1));
  }
  return res;
}
