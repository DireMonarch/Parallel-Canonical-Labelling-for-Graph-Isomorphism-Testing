# Copyright 2025 Jim Haslett

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import networkx as nx
# import numpy as np

__all__ = ['johnson']

@nx._dispatchable(graphs=None, returns_graph=True)
def johnson_graph(n, k):
    """Returns the adjacency matrix of a johnson graph J(n,k).
       For: n > 0, k > 0, k < n."""
    if n < 1 or k < 1 or k > n:
        return None
    subsets = list(filter(lambda x: len(x) == k, _powerset(range(n))))

    size = len(subsets)
    G : nx.graph.Graph = nx.empty_graph(size)

    # adjacency = np.zeros((size, size))
    for i in range(size):
        for j in range(i+1, size):
            if len(subsets[i].intersection(subsets[j])) == k-1:
                G.add_edge(i,j)
                
    #             adjacency[i][j] = 1
    #             adjacency[j][i] = 1
    return G
    # return adjacency


def _powerset(A):
    """Returns the powerset of a given set."""
    P = [set([])]
    for i in range(len(A)):
        for j in range(len(P)):
            subset = P[j].copy()
            subset.add(i)
            P.append(set(subset))
    return P
