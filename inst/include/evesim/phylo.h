//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT


#ifndef PHYLO_H_INCLUDED
#define PHYLO_H_INCLUDED
#pragma once

#include <climits>  // CHAR_BIT
#include "evesim/rutils.h"


namespace tres_sim { 

  // ape::phylo wrapper 
  class phylo_t {
  public:
    explicit phylo_t(Rcpp::List&& phy);
    explicit phylo_t(const struct tree_t& tree);

    auto e0() const noexcept { return edge_.column(0); }
    auto e1() const noexcept { return edge_.column(1); }
    auto len() const noexcept { return edge_length_; }
    auto tip_label() const noexcept { return tip_label_; }

    auto e0() noexcept { return edge_.column(0); }
    auto e1() noexcept { return edge_.column(1); }
    auto len() noexcept { return edge_length_; }
    auto tip_label() noexcept { return tip_label_; }

    int size() const noexcept { return len().size(); }
    int tips() const noexcept { return tip_label_.size(); }
    int nnode() const noexcept { return tips() - 1; }

    // moves ape::phylo object out, destroys internals
    Rcpp::List unwrap();

  private:
    Rcpp::List phy_;
    RcppParallel::RMatrix<int> edge_;
    RcppParallel::RVector<double> edge_length_;
    RcppParallel::RVector<int> tip_label_;
  };


} // namespace tres_sim


#endif
