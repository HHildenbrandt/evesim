//  SPDX-FileCopyrightText: 2024 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
//  SPDX-License-Identifier: MIT

#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED
#pragma once

#include <cmath>
#include <cassert>
#include <vector>
#include <array>
#include <numeric>


namespace tres_sim {

  // used to rectify extremal times in transformed trees.
  constexpr auto present_tol = 7.0 / (1000000.0 * 364.0);   // one week in 1ma


  // minimal node type for tree representation.
  struct node_t {
    double t;                 // tend
    int ances;                // index into node table, -1 for root
    std::array<int, 2> desc;  // indices into node table, {-1, -1} for leafs
    int label;                // carries over the signed index from ltable_t, one-based; undefined for inodes

    bool is_root() const noexcept { return ances == -1; }
    bool is_inode() const noexcept { return (desc[0] >= 0); }         // binary tree
    bool is_leaf() const noexcept { return (desc[0] < 0); }           // binary tree
    bool is_deg() const noexcept { return (desc[0] * desc[1]) < 0; }  // degenerated node

    // solves for 'i' in node.desc[i] == idx.
    int desc_idx(int idx) const noexcept { return static_cast<int>(desc[0] != idx); }
    
    friend bool operator==(const node_t& a, const node_t& b) noexcept {    // for debugging
      return (a.t == b.t)
          && (a.ances == b.ances)
          && (a.desc == b.desc)
          && (a.label == b.label);
    }
  };


  // tree representation.
  struct tree_t {
    double age = 0.0;
    int root = 0;
    bool ultrametric = false;    // extinct species dropped
    std::vector<node_t> nodes;

    static tree_t from(const class ltable_t&, bool drop_extinct);
    static tree_t from(const class phylo_t&);

    int tips() const noexcept { return root; }
    int nnode() const noexcept { return static_cast<int>(nodes.size()) - root; }
    void clear() noexcept { age = 0.0; root = 0; ultrametric = false; nodes.clear(); }

    template <typename OIT>
    void tip_label(OIT oit) const {
      for (int i = 0; i < tips(); ++i) {
        *oit++ = std::abs(nodes[i].label);
      }
    }
  };


  namespace detail {

    // helper for drop_tips_if: skip over nodes[idx].ances
    inline void drop_tip(std::vector<node_t>& nodes, int idx) {
      const auto ances = nodes[idx].ances;
      auto& node_ances = nodes[ances];
      const auto aa = node_ances.ances;
      assert(aa != -1);
      auto& node_aa = nodes[aa];
      const auto adesc = node_aa.desc_idx(ances);
      node_aa.desc[adesc] = node_ances.desc[!node_ances.desc_idx(idx)];
      nodes[node_aa.desc[adesc]].ances = aa;
    }


    // helper for drop_tips_if: reorder by node_map
    // return true if ultrametric
    inline bool reorder(std::vector<node_t>& nodes, int src_tips, std::vector<int>& node_map) {
      int next_node = 0;
      for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (node_map[i] == 0) {
          node_map[i] = next_node++;
        }
      }
      bool ultrametric = true;
      for (int i = 0; i < src_tips; ++i) {
        if (node_map[i] != -1) {
          const auto& node = nodes[i];
          nodes[node_map[i]] = { node.t, node_map[node.ances], { -1,-1 }, node.label };
          ultrametric &= (node.t == 0.0);
        }
      }
      for (int i = src_tips; i < static_cast<int>(nodes.size()); ++i) {
        if (node_map[i] != -1) {
          const auto& node = nodes[i];
          nodes[node_map[i]] = { node.t, node_map[node.ances], { node_map[node.desc[0]], node_map[node.desc[1]] }, 0 };
        }
      }
      nodes.resize(next_node);
      return ultrametric;
    }

  }


  template <typename Pred>
  inline void drop_tips_if_inplace(tree_t& tree, Pred pred) {
    auto& nodes = tree.nodes;
    auto node_map = std::vector<int>(nodes.size(), 0);
    for (int i = 0; i < tree.tips(); ++i) {
      auto& node = nodes[i];
      if (pred(node)) {
        // ancestor has become degenerated, skip over it
        detail::drop_tip(nodes, i);
        node_map[i] = node_map[node.ances] = -1;   // mark as invalid
      }
    }
    tree.ultrametric = detail::reorder(nodes, tree.tips(), node_map);
    tree.root = node_map[tree.root];
    nodes[tree.root].ances = -1;
  }


  template <typename Pred>
  inline tree_t drop_tips_if(const tree_t& src, Pred&& pred) {
    auto tree = src;
    drop_tips_if_inplace(tree, std::forward<Pred>(pred));
    return tree;
  }

}

#endif