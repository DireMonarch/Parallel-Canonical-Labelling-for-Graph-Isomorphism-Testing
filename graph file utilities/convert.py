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
import os
import networkx as nx



def read_dimacs(filename: str) -> nx.graph.Graph:
    G = nx.Graph()
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('c'):
                continue
            if line.startswith('p'):
                parts = line.split()
                if len(parts) != 4 or parts[1] != 'edge':
                    raise ValueError('Invalid DIMACS format')
                n = int(parts[2])
                m = int(parts[3])
                G.add_nodes_from(range(1, n+1))
            elif line.startswith('e'):
                parts = line.split()
                if len(parts) != 3:
                    raise ValueError('Invalid DIMACS format')
                u = int(parts[1])
                v = int(parts[2])
                G.add_edge(u, v)
    return G


def usage():
    print('Usage: ')
    print(f'{sys.argv[0]} -e <filename>  \n\tConverts an edge list in <filename> to sparse6 format')
    print(f'{sys.argv[0]} -d <filename>  \n\tConverts a DIMACS file in <filename> to sparse6 format')


def main():
    if len(sys.argv) < 3:
        usage()
        exit(0)

    mode = sys.argv[1]
    if mode not in ['-e', '-d']:
        usage()
        exit(0)
        
    infilename = sys.argv[2]
    
    outfilename, ext = os.path.splitext(os.path.basename(infilename))
    if ext in ['.gz', '.bz2', '.xz']:
        outfilename = os.path.splitext(outfilename)[0]
    outfilename += '.s6'
    print(outfilename)
    
    if mode == '-e':
        G : nx.graph.Graph = nx.read_edgelist(infilename)
        nx.sparse6.write_sparse6(G, outfilename)
        print(f'Read edge list from {infilename} and wrote sparse6 format to {outfilename}')


    if mode == '-d':
        G : nx.graph.Graph = read_dimacs(infilename)
        nx.sparse6.write_sparse6(G, outfilename)
        print(f'Read DIMACS from {infilename} and wrote sparse6 format to {outfilename}')


if __name__ == '__main__':
    main()    