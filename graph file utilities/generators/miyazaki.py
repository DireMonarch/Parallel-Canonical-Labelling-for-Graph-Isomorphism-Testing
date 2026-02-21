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

__all__ = ['miyazaki']

@nx._dispatchable(graphs=None, returns_graph=True)
def miyazaki(k):
    if k < 1:
        raise nx.NetworkXError("must satisfy k >= 1")
    
    n = k * 20
    G : nx.graph.Graph = nx.empty_graph(n)
    
    #
    # Always add from smaller vertex number to larger to avoid duplicate edges
    #
    G.add_edge(0,3)
    G.add_edge(0,4)
    G.add_edge(0,6)
    
    G.add_edge(1,2)
    G.add_edge(1,5)
    G.add_edge(1,7)
    
    G.add_edge(2,5)
    G.add_edge(2,6)
    
    G.add_edge(3,4)
    G.add_edge(3,7)
    
    G.add_edge(4,8)
    G.add_edge(5,8)
    
    G.add_edge(6,9)
    G.add_edge(7,9)
    
    #
    # 10 and 11 are in the next "section"
    # Always bridge between sections from smaller to larter vertex number
    #   to avoid duplicate eges
    #
    G.add_edge(8,10)
    G.add_edge(9,11)    
    
    
    for i in range(k-1):
        b = (i * 20) + 10  # base (zero) vertex for this bridge section
        
        #
        # Left half of bridge section
        #
        G.add_edge(b+0, b+2)    # 10, 12
        G.add_edge(b+0, b+3)    # 10, 13
        
        G.add_edge(b+1, b+4)    # 11, 14
        G.add_edge(b+1, b+5)    # 11, 15
        
        G.add_edge(b+2, b+6)    # 12, 16
        G.add_edge(b+2, b+8)    # 12, 18
        
        G.add_edge(b+3, b+7)    # 13, 17
        G.add_edge(b+3, b+9)    # 13, 19
        
        G.add_edge(b+4, b+7)    # 14, 17
        G.add_edge(b+4, b+8)    # 14, 18
        
        G.add_edge(b+5, b+6)    # 15, 16
        G.add_edge(b+5, b+9)    # 15, 19
        
        #
        # Bridge section right to left connectors
        #
        G.add_edge(b+6, b+10)   # 16, 20
        G.add_edge(b+7, b+11)   # 17, 21
        G.add_edge(b+8, b+12)   # 18, 22
        G.add_edge(b+9, b+13)   # 19, 23
        
        #
        # Right half of bridge section
        #
        G.add_edge(b+10, b+14)  # 20, 24
        G.add_edge(b+10, b+16)  # 20, 26
        
        G.add_edge(b+11, b+15)  # 21, 25
        G.add_edge(b+11, b+17)  # 21, 27
        
        G.add_edge(b+12, b+15)  # 22, 25
        G.add_edge(b+12, b+16)  # 22, 26
        
        G.add_edge(b+13, b+14)  # 23, 24
        G.add_edge(b+13, b+17)  # 23, 27
        
        G.add_edge(b+14, b+18)  # 24, 28
        G.add_edge(b+15, b+18)  # 25, 28
        
        G.add_edge(b+16, b+19)  # 26, 29
        G.add_edge(b+17, b+19)  # 27, 29
        
        #
        # link to next section
        #
        G.add_edge(b+18, b+20)  # 28, 30
        G.add_edge(b+19, b+21)  # 28, 31
        
    
    
    #
    # add Rigth circuit
    #
    b = ((k-1) * 20) + 10  # base (zero) vertex for this bridge section
    #
    # Left half of bridge section
    #
    G.add_edge(b+0, b+2)    # 10, 12
    G.add_edge(b+0, b+3)    # 10, 13
    
    G.add_edge(b+1, b+4)    # 11, 14
    G.add_edge(b+1, b+5)    # 11, 15
    
    G.add_edge(b+2, b+6)    # 12, 16
    G.add_edge(b+2, b+8)    # 12, 18
    
    G.add_edge(b+3, b+7)    # 13, 17
    G.add_edge(b+3, b+9)    # 13, 19
    
    G.add_edge(b+4, b+7)    # 14, 17
    G.add_edge(b+4, b+8)    # 14, 18
    
    G.add_edge(b+5, b+6)    # 15, 16
    G.add_edge(b+5, b+9)    # 15, 19
    
    G.add_edge(b+6, b+9)    # 16, 19
    G.add_edge(b+7, b+8)    # 17, 18
    
    return G

