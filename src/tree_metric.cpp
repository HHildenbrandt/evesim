//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#include <tuple>
#include <stack>
#include <deque>
#include <numeric>
#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include "evesim/tree_metric.h"
#include "evesim/rutils.h"


namespace tres_sim {

  namespace tree_metric {

    void set_dim_names(const Rcpp::RObject& MV, const tree_t& tree) {
      auto names = Rcpp::IntegerVector(tree.tips());
      tree.tip_label(names.begin());
      if (Rcpp::is<Rcpp::NumericMatrix>(MV)) {
        Rcpp::as<Rcpp::NumericMatrix>(MV).attr("dimnames") = Rcpp::List::create(names, names);
        return;
      }
      if (Rcpp::is<Rcpp::NumericVector>(MV)) {
        Rcpp::as<Rcpp::NumericVector>(MV).attr("names") = names;
        return;
      }
      // silent fail
    }


    cophenetic::dij_t::dij_t(const tree_t& tree) : tree(tree), visitors(tree.nnode(), -1) {
      const auto& nodes = tree.nodes;
      auto* first_visitor = visitors.data() - tree.root; // iterator would be UB
      // first_visitor[n] == x <-> tip x is the first one which path to root traverses node n 
      for (int i = 0; i < tree.tips(); ++i) {
        auto ai = nodes[i].ances;
        while ((ai != -1) && (first_visitor[ai] == -1)) {
          first_visitor[ai] = i;
          ai = nodes[ai].ances;
        }
      }
    }


    // inner loop body for `apply`
    // returns distance between tips i and j
    //
    // this is not a standalone solution because the outer loop
    // needs to match a specific sequencing.
    double cophenetic::dij_t::operator()(double const* const Di, int i, int j) const {
      const auto& nodes = tree.nodes;
      if ((nodes[i].label > 0) != (nodes[j].label > 0)) {
        // different root clade
        return 2.0 * nodes[tree.root].t - (nodes[i].t + nodes[j].t);
      }
      int ai = nodes[i].ances;
      int aj = nodes[j].ances;
      const auto* first_visitor = visitors.data() - tree.root; // iterator would be UB
      while (ai != aj) {
        if (nodes[ai].t > nodes[aj].t) {
          if (auto other = first_visitor[aj]; other != j) {
            // Di[other] seen in earlier iteration
            return Di[other] + ((nodes[other].t) - (nodes[j].t));
          }
          aj = nodes[aj].ances;
        }
        else {
          ai = nodes[ai].ances;
        }
      }
      // ai is mra
      return 2.0 * nodes[ai].t - (nodes[i].t + nodes[j].t);
    }


    // pairwise distance between tips, generic storage
    void cophenetic::apply(const tree_t& tree, RcppParallel::RMatrix<double> D) {
      rutils::tbb_global_control_guard gc(false);
      auto dij = dij_t(tree);
      if (gc.num_threads() <= 2) {
        // serial
        for (int i = 0; i < tree.tips(); ++i) {
          double* Di = D.begin() + i * tree.tips();
          Di[i] = 0.0;
          for (int j = i + 1; j < tree.tips(); ++j) {
            D(i, j) = D(j, i) = dij(Di, i, j);
          }
        }
      }
      else {
        // parallel, more work
        tbb::parallel_for(tbb::blocked_range<int>(0, tree.tips()), [&](const auto& r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            double* Di = D.begin() + i * tree.tips();
            for (int j = 0; j < i; ++j) {
              Di[j] = dij(Di, i, j);
            }
            Di[i] = 0.0;
            for (int j = i + 1; j < tree.tips(); ++j) {
              Di[j] = dij(Di, i, j);
            }
          }
        });
      }
    }


    // pairwise distance between tips
    Rcpp::NumericMatrix cophenetic::operator()(const tree_t& tree) const {
      auto DR = Rcpp::NumericMatrix(tree.root, tree.root);
      apply(tree, RcppParallel::RMatrix<double>(DR));
      return DR;
    }


    // mean evolutionary distinctiveness, generic storage
    void ed::apply(const tree_t& tree, RcppParallel::RVector<double> D) {
      const int tips = tree.tips();
      auto UDD = std::unique_ptr<double[]>(new double[tips * tips]);
      auto dd = rutils::RView(rutils::RView(UDD.get(), tips, tips));
      cophenetic::apply(tree, dd);
      rutils::tbb_global_control_guard gc{false};
      if (gc.num_threads() == 1) {
        for (int i = 0; i < tips; ++i) {
          const auto* first = (dd.begin() + i * tips);
          D[i] = std::accumulate(first, first + tips, 0.0) / (tips - 1);
        }
      }
      else {
        tbb::parallel_for(tbb::blocked_range<int>(0, tips), [&](const auto& r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            const auto* first = (dd.begin() + i * tips);
            D[i] = std::accumulate(first, first + tips, 0.0) / (tips - 1);
          }
        });
      }
    }


    // mean evolutionary distinctiveness
    Rcpp::NumericVector ed::operator()(const tree_t& tree) const {
      auto D = Rcpp::NumericVector(tree.root);
      apply(tree, RcppParallel::RVector<double>(D));
      return D;
    }


    // nearest neighbor distance between tips, generic storage
    void nnd::apply(const tree_t& tree, RcppParallel::RVector<double> D) {
      const auto& nodes = tree.nodes;
      const int tips = tree.tips();
      if (tree.ultrametric) {
        // all leafs are at age == 0.0...
        // not enough work for parallel execution
        for (int i = 0; i < tips; ++i) {
          const auto a = nodes[i].ances;
          D[i] = 2.0 * nodes[a].t;
        }
        return;
      }
      auto UDD = std::unique_ptr<double[]>(new double[tips * tips]);
      auto dd = rutils::RView(rutils::RView(UDD.get(), tips, tips));
      cophenetic::apply(tree, dd);
      rutils::tbb_global_control_guard gc{};
      tbb::parallel_for(tbb::blocked_range<int>(0, tips), [&](const auto& r) {
        for (int i = r.begin(); i != r.end(); ++i) {
          double* Di = dd.begin() + i * tips;
          Di[i] = std::numeric_limits<double>::max();
          D[i] = *std::min_element(Di, Di + tips);
        }
      });
    }
    

    Rcpp::NumericVector nnd::operator()(const tree_t& tree) const {
      auto DR = Rcpp::NumericVector(tree.root);
      apply(tree, RcppParallel::RVector<double>(DR));
      return DR;
    }


    void pd::apply(const tree_t& tree, double& D) {
      D = 0.0;
      const auto& nodes = tree.nodes;
      for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (i != tree.root) {
          auto& node = nodes[i];
          D += nodes[node.ances].t - node.t;
        }
      }
    }


    double pd::operator()(const tree_t& tree) const {
      double D;
      apply(tree, D);
      return D;
    }


    void mpd::apply(const tree_t& tree, double& D) {
      const auto tips = tree.tips();
      auto UDD = std::unique_ptr<double[]>(new double[tips * tips]);
      auto dd = rutils::RView<double>(UDD.get(), tips, tips);
      cophenetic::apply(tree, dd);

      rutils::tbb_global_control_guard gc{false};
      double sum = 0.0;
      if (gc.num_threads() == 1) {
        sum = std::accumulate(dd.begin(), dd.end(), 0.0);
      }
      else {
        sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, dd.size()), 0.0, [&](const auto& r, double running_sum) {
          running_sum += std::accumulate(dd.begin() + r.begin(), dd.begin() + r.end(), 0.0);
          return running_sum;
        }, std::plus<double>() );
      }
      // lower.tri
      D = sum / (tips * (tips - 1));
    }


    double mpd::operator()(const tree_t& tree) const {
      double D;
      apply(tree, D);
      return D;
    }

  };  // namespace tree_metric

} // namespace tres_sim
