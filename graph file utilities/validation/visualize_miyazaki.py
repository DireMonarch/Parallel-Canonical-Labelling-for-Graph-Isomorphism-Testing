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



def generate_positions(k: int):
    pos = {0:(-.25,8), 1:(0,6), 2:(0,2), 3:(-.25,0),
        4:(1,7), 5:(1,5), 6:(1,3), 7:(1,1),
                    8:(2,5), 9:(2,3)}
    
    
    p = 10
    x = 2
    for i in range(1, k):
        #         10:(3,5), 11:(3,3),
        x += 1
        for y in [5, 3]:
            pos[p] = (x, y)
            p += 1
        
        # 12:(4,7), 13:(4,5), 14:(4,3), 15:(4,1),
        x += 1
        for y in [7, 5, 3, 1]:
            pos[p] = (x, y)
            p += 1
                    
        # 16:(5,8), 17:(5,6), 18:(5,2), 19:(5,0),
        x += 1
        for y in [8, 6, 2, 0]:
            pos[p] = (x, y)
            p += 1
            
        # 20:(6,8), 21:(6,6), 22:(6,2), 23:(6,0),
        x += 1
        for y in [8, 6, 2, 0]:
            pos[p] = (x, y)
            p += 1
            
        # 24:(7,7), 25:(7,5), 26:(7,3), 27:(7,1),
        x += 1
        for y in [7, 5, 3, 1]:
            pos[p] = (x, y)
            p += 1
            
        #             28:(8,5), 29:(8,3),
        x += 1
        for y in [5, 3]:
            pos[p] = (x, y)
            p += 1
        
        
        
    
    #         30:(9,5), 31:(9,3),
    x += 1
    for y in [5, 3]:
        pos[p] = (x, y)
        p += 1
        
    # 32:(10,7), 33:(10,5), 34:(10,3), 35:(10,1),
    x += 1
    for y in [7, 5, 3, 1]:
        pos[p] = (x, y)
        p += 1
        
    # 36:(11.25,8), 37:(11,6), 38:(11,2), 39:(11.25,0)}
    x += 1
    for y in [8, 6, 2, 0]:
        if y in [8, 0]:
            pos[p] = (x+.25, y)
        else:
            pos[p] = (x, y)
        p += 1

    
    
    return pos


def visualize(G: nx.graph.Graph, k: int):
    pos = generate_positions(k)
    plt.figure(figsize=(6.4*k, 4.8))
    nx.draw_networkx(G, pos=pos)
    plt.savefig(f'M-{k} generated.png')
    # plt.show()

def usage():
    print('Usage: ')
    print(f'{sys.argv[0]} <k> <filename>  \n\tGenerates Miyazaki graph of size k | k > 0')


def main():
    
    if len(sys.argv) < 3:
        usage()
        exit(0)

    k = int(sys.argv[1])
    if k < 1:
        usage()
        exit(1)

    filename = sys.argv[2]

    visualize(nx.read_graph6(filename), k)


if __name__ == '__main__':
    main()  