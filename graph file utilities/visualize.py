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
import matplotlib.pyplot as plt




def usage():
    print('Usage: ')
    print(f'{sys.argv[0]} <filename>')
    print('\tReads a graph in sparse6 format from <filename> and generates a PNG image of the graph')


def main():
    
    if len(sys.argv) < 2:
        usage()
        exit(0)


    filename = sys.argv[1]
    
    G = nx.read_sparse6(filename)
    plt.figure(figsize=(6.4, 4.8))
    nx.draw_networkx(G)
    plt.savefig(f'{filename}.png')
    # plt.show()

if __name__ == '__main__':
    main()  