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

import sys
import networkx as nx
from generators.miyazaki import miyazaki
from generators.johnson import johnson_graph

def usage():
    print('Usage: ')
    print(f'{sys.argv[0]} -M <k> <filename>  \n\tGenerates Miyazaki graph of size k')
    print(f'{sys.argv[0]} -J <n> <k> <filename>  \n\tGenerates Johnson graph of type J(n, k)')


def main():
    if len(sys.argv) < 4:
        usage()
        exit(0)

    if sys.argv[1] == '-M':
        k = int(sys.argv[2])
        filename = sys.argv[3]
        
        G : nx.graph.Graph = miyazaki(k)
        print(f'Generated Miyazaki Graph with K={k}:  {G}')
        nx.graph6.write_graph6(G, filename)
        print(f'Written to filename: {filename}')
        
    elif sys.argv[1] == "-J" and len(sys.argv) > 4:
        n = int(sys.argv[2])
        k = int(sys.argv[3])
        filename = sys.argv[4]
                
        
        G : nx.graph.Graph = johnson_graph(n,k)
        
        print(f'Generated Johnson Graph of type J({n}, {k}):  {G}')
        nx.graph6.write_graph6(G, filename)
        print(f'Written to filename: {filename}')
        
    else:
        usage()

if __name__ == '__main__':
    main()    