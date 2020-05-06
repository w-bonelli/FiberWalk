"""
FIBER WALK DEMO

The Fiber Walk class to run a Fiber Walk simulation.

The demo uses:
- the networkX package (http://networkx.lanl.gov/)
- the numpy package (http://sourceforge.net/projects/numpy/)
- the scipy package (http://www.scipy.org/SciPy)
- the classes Lattice and FiberWalk written by Alexander Bucksch

Note: The graphics shown in the paper where generated with mayavi2
(http://docs.enthought.com/mayavi/mayavi/index.html)

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the Fiber Walk Paper if you use the code for your scientific project.

-------------------------------------------------------------------------------------------
Author: Alexander Bucksch
School of Biology and Interactive computing
Georgia Institute of Technology

Mail: bucksch@gatech.edu
Web: http://www.bucksch.nl
-------------------------------------------------------------------------------------------

Copyright (c) 2012 Alexander Bucksch
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the Fiber Walk Demo Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
# !/usr/bin/python

import lattice as L
import os
import pickle
import numpy as np
import networkx as nx


class FiberWalk():

    def __init__(self, dim, objects=1, contractions=1, min_length=0):
        self.__lattice = L.Lattice(dim)  # the lattice
        self.__walk = nx.Graph()  # The walk itself as a networkX object
        self.__dim = dim  # dimension of the walk
        self.__min_length = min_length
        self.__contractions = contractions
        self.__positions = []
        self.__distances_to_origin = []
        self.__all_distances_to_origin = []
        self.__valid_neighbors = []
        self.__contract_neighbors = []
        self.__sa = []
        self.__sa_path = []
        self.__contr = []
        self.__path_labels = []
        self.__objects = objects  # number of walks to be computed
        self.__all_walk_distances = []  # distances to the origin for all walks computed
        self.__all_valid_neighbors = []
        self.__all_contract_neighbors = []
        self.__all_sa = []
        self.__all_sa_paths = []
        self.__all_contr = []
        self.__all_walks_sa = []
        self.__all_walks_sa_paths = []
        self.__allWalksContr = []
        self.__norm = 2  # norm used to calculate the distance from the seed
        self.__sa_count = 0
        self.__sa_path_count = 0
        self.__contr_count = 0

        # check for directories to enable restart
        try:
            os.stat('./tmp/')
        except:
            os.makedirs('./tmp/')

    # getters

    def get_lattice(self):
        return self.__lattice

    def get_walk(self):
        return self.__walk

    def get_all_distances_to_origin(self):
        return self.__all_distances_to_origin

    # reset functions

    def reset_arrays(self, reset_lattice=False):
        if reset_lattice:
            self.__lattice.reset_Lattice()
        self.__walk.clear()
        self.__positions = []
        self.__distances_to_origin = []
        self.__valid_neighbors = []
        self.__contract_neighbors = []
        self.__path_labels = []

    def reset_all_arrays(self):
        self.__all_distances_to_origin = []
        self.__all_valid_neighbors = []
        self.__all_contract_neighbors = []
        self.__all_sa = []
        self.__all_sa_paths = []
        self.__all_contr = []

    def define_seed(self, pos='center'):
        ret = np.zeros(self.__dim)
        if self.__lattice != 0:
            n = self.__lattice.get_Lattice().node[tuple(list(ret))]
            # mark the seed as consumed
            n['consumed'] = True
            n['counter'] = 0
            n['pos'] = ret
            n['SA'] = False
            n['Contr'] = 1
        return tuple(list(ret))

    def select_walk_edge_leading_label_self_avoid_contract_sq(self, new_pos, old_pos, new_walk):
        stop = False
        neighbors = self.__lattice.get_Lattice().neighbors(new_pos)
        back_edge = self.__lattice.get_Lattice()[new_pos][old_pos][0]['vdir']
        # COLLECT NEIGHBOURS AND EXCLUDE BACK EDGES AND SELFAVOIDING EDGES
        valid_neighbors = []
        for i in neighbors:
            if not self.__lattice.get_Lattice()[new_pos][i][0]['SA']:
                if self.__walk.has_edge(new_pos, i):
                    print 'BUG !!!!!!'
                e = self.__lattice.get_Lattice()[new_pos][i][0]['vdir']
                # check if the walk does not form a loop -> self avoidance
                if not self.__lattice.get_Lattice().node[i]['consumed']:
                    test_back = np.zeros_like(back_edge)
                    for k in range(len(back_edge)):
                        test_back[k] = back_edge[k] + e[k]
                        test_back[k] = np.abs(test_back[k])
                    test_back = np.sum(test_back)
                    if test_back <= self.__dim:
                        valid_neighbors.append(i)
                else:
                    if not self.__lattice.get_Lattice()[new_pos][i][0]['SA']:
                        if not self.__walk.has_edge(new_pos, i):
                            self.__sa_path_count += 1
                        self.__lattice.get_Lattice()[new_pos][i][0]['SA'] = True
                        self.__lattice.get_Lattice()[i][new_pos][0]['SA'] = True

            elif not new_walk:
                if self.__walk.has_edge(i, new_pos):
                    print 'BUG !!!!!!'

        if len(valid_neighbors) == 0:
            return new_pos, new_pos, True, False

        old_pos = new_pos

        r_idx = np.random.randint(0, len(valid_neighbors));
        new_pos = valid_neighbors[r_idx]

        # check for valid contractions -> Just Degug code
        for i in range(self.__contractions):
            self.get_lattice().expandLattice(new_pos)
            contr_ok = self.contract(new_pos, old_pos)
            if not contr_ok:
                return new_pos, new_pos, True, False

        self.get_lattice().expandLattice(new_pos)
        if len(valid_neighbors) != 0:
            e = self.__lattice.get_Lattice()[old_pos][new_pos][0]['vdir']  # get the path label
            self.__path_labels.append(e)  # store the path labels
            # add edge to the walk
            if not self.__walk.has_edge(old_pos, new_pos):
                self.__walk.add_node(new_pos, counter=1, pos=None, trapped=False)
                self.__walk.add_edge(new_pos, old_pos, counter=1)
                self.__lattice.get_Lattice()[new_pos][old_pos][0]['SA'] = True
                self.__lattice.get_Lattice()[old_pos][new_pos][0]['SA'] = True
                self.__lattice.get_Lattice()[new_pos][old_pos][0]['consumed'] = True
                self.__lattice.get_Lattice()[old_pos][new_pos][0]['consumed'] = True
                new_walk = True
            else:
                self.__walk[old_pos][new_pos]['counter'] += 1
                self.__walk.node[new_pos]['counter'] += 1
                new_walk = False

            self.__valid_neighbors.append(len(valid_neighbors))
            # Set vertex as consumed
            self.__lattice.get_Lattice().node[new_pos]['consumed'] = True;
        else:
            stop = True

        return new_pos, old_pos, stop, new_walk

    def init_seed_contraction(self, newPos, oldPos):
        # Contract the Seed
        new_walk = False
        valid_neighbors = []
        path_idx = []
        contraction_neighbors = []
        # Consume the seed ->Don't consume it, otherwise it is not possible to walk out there
        self.__lattice.get_Lattice().node[oldPos]['consumed'] = True
        # contract the edges of the newly reached vertex
        neighbors = self.__lattice.get_Lattice().neighbors(tuple(oldPos))
        for i in neighbors:
            if self.__walk.has_edge(oldPos, i):
                path_idx.append(i)
            elif not self.__lattice.get_Lattice().node[i]['consumed']:
                valid_neighbors.append(i)
            elif self.__lattice.get_Lattice().node[i]['consumed']:
                if not self.__walk.has_edge(oldPos, i):
                    self.__sa_path_count += 1
                self.__lattice.get_Lattice()[newPos][i][0]['SA'] = True
                self.__lattice.get_Lattice()[i][newPos][0]['SA'] = True

        if len(valid_neighbors) > 0:
            rIdx = np.random.randint(0, len(valid_neighbors));
            newPos = valid_neighbors[rIdx]

            # make array of S to increment
            newPos = np.array(newPos)
            newPos = tuple(list(newPos))
            if not self.__lattice.get_Lattice().node[newPos]['consumed']:
                new_walk = True
            else:
                new_walk = False
            self.get_lattice().expandLattice(newPos)
            e = np.array(oldPos) - np.array(newPos)
            neighbors = self.__lattice.get_Lattice().neighbors(newPos)

            self.__lattice.get_Lattice().node[newPos]['consumed'] = True

            if not self.__walk.has_edge(oldPos, newPos):
                self.__walk.add_node(newPos, counter=1, pos=None, trapped=False)
                self.__walk.add_edge(oldPos, newPos, counter=1)
                self.__lattice.get_Lattice()[newPos][oldPos][0]['SA'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['SA'] = True
                self.__lattice.get_Lattice()[newPos][oldPos][0]['consumed'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['consumed'] = True
                new_walk = True
            else:
                self.__walk[oldPos][newPos]['counter'] += 1
                self.__walk.node[newPos]['counter'] += 1
                new_walk = False

            self.__path_labels.append(e)
            self.__contract_neighbors.append(len(contraction_neighbors))
            # all neighbors are valid in the init procedure
            self.__valid_neighbors.append(len(list(neighbors)))
            stop = False
        else:
            stop = True

        return newPos, oldPos, stop, new_walk

    # debug function to test if a contraction was valid
    def contraction_is_valid(self, new_pos, i, j):
        ret = True

        # trivial case
        if j == new_pos:
            ret = False
        elif i == new_pos:
            ret = False
        elif i == j:
            ret = False
        elif self.__lattice.get_Lattice().node[i]['consumed']:
            if not self.__walk.has_edge(new_pos, i):
                self.__lattice.get_Lattice()[new_pos][i][0]['SA'] = True
                self.__lattice.get_Lattice()[i][new_pos][0]['SA'] = True
                ret = False
        # detect a possible triangle
        elif self.__lattice.get_Lattice().has_edge(new_pos, j):
            if self.__lattice.get_Lattice().has_edge(new_pos, i):
                l1 = self.__lattice.get_Lattice()[new_pos][i][0]['vdir']
                l2 = self.__lattice.get_Lattice()[new_pos][j][0]['vdir']
                # check label equality
                test_eq = True
                for k in range(len(l1)):
                    if l1[k] != l2[k]:
                        test_eq = False
                        break
                if test_eq:
                    ret = False
                    if not self.__lattice.get_Lattice()[new_pos][i][0]['SA']:

                        self.__lattice.get_Lattice()[new_pos][i][0]['SA'] = True
                        self.__lattice.get_Lattice()[i][new_pos][0]['SA'] = True
                        if not self.__walk.has_edge(new_pos, i): self.__sa_count += 1
                    if not self.__lattice.get_Lattice()[new_pos][j][0]['SA']:

                        self.__lattice.get_Lattice()[new_pos][j][0]['SA'] = True
                        self.__lattice.get_Lattice()[j][new_pos][0]['SA'] = True
                        if not self.__walk.has_edge(new_pos, j):  self.__sa_count += 1

                else:
                    ret = True

        if ret:
            if self.__lattice.get_Lattice().has_edge(i, j):
                if self.__lattice.get_Lattice().has_edge(i, new_pos):
                    l1 = self.__lattice.get_Lattice()[i][new_pos][0]['vdir']
                    l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                    test_eq = True
                    for k in range(len(l1)):
                        if l1[k] != l2[k]:
                            test_eq = False
                            break
                    if test_eq:
                        ret = False
                        if not self.__lattice.get_Lattice()[j][i][0]['SA']:
                            self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                            self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                            self.__sa_count += 1
                        if not self.__lattice.get_Lattice()[i][j][0]['SA']:
                            self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                            self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                            self.__sa_count += 1
            else:
                ret = True
        if ret:
            if self.__lattice.get_Lattice().has_edge(j, new_pos):
                l1 = self.__lattice.get_Lattice()[new_pos][j][0]['vdir']
                l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                test_eq = True
                for k in range(len(l1)):
                    if l1[k] != l2[k]:
                        test_eq = False
                        break
                if test_eq == True:
                    ret = False
                    if not self.__lattice.get_Lattice()[j][i][0]['SA']:
                        self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                        self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                        if not self.__walk.has_edge(new_pos, i): self.__sa_count += 1
                    if not self.__lattice.get_Lattice()[j][i][0]['SA']:
                        self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                        self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                        if not self.__walk.has_edge(new_pos, i): self.__sa_count += 1
        else:
            ret = True
        return ret

    # funtion to perform the contraction     
    def contract(self, new_pos, old_pos):
        # Start Contraction
        contraction_neighbors = []
        # get the neighbors of the random selected vertex
        neighbors = self.__lattice.get_Lattice().neighbors(new_pos)
        e = self.__lattice.get_Lattice()[old_pos][new_pos][0]['vdir']
        # collect all neighbors that are possible to contract
        for i in neighbors:
            e2 = self.__lattice.get_Lattice()[new_pos][i][0]['vdir']
            testEq = True
            #            check the back edge of the walk
            for n in range(len(e)):
                if e[n] != e2[n]:
                    testEq = False
            #            if it is not the back edge
            if not testEq:
                # check for self avoidance
                if not self.__lattice.get_Lattice().node[i]['consumed']:
                    contraction_neighbors.append(i)
                    self.__contr_count += 1
                # check if there is space otherwise move the walk to the opposite position
                elif self.__lattice.get_Lattice()[new_pos][i][0]['Contr'] < self.__lattice.get_Lattice().node[new_pos][
                    'counter']:
                    return False

        for i in contraction_neighbors:
            n2 = list(self.__lattice.get_Lattice().neighbors(i))
            l1 = self.__lattice.get_Lattice()[new_pos][i][0]['vdir']
            self.get_lattice().expandLattice(i)

            for j in n2:
                if j != new_pos:
                    l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                    l = l1 + l2
                    for a in range(len(l)):
                        if l[a] == 0:
                            l[a] = 0
                        if l[a] < 0:
                            l[a] = -1
                        if l[a] > 0:
                            l[a] = 1

                    if self.contraction_is_valid(new_pos, i, j):

                        # Note: This will make every edge of the walk a non contractable edge
                        c = self.__lattice.get_Lattice()[i][j][0]['Contr']
                        if not self.__lattice.get_Lattice().has_edge(new_pos, j):
                            self.__lattice.get_Lattice().add_edge((new_pos), j, vdir=l, SA=False, Contr=c + 1)
                        if not self.__lattice.get_Lattice().has_edge(j, new_pos):
                            self.__lattice.get_Lattice().add_edge(j, (new_pos), vdir=(l * (-1)), SA=False, Contr=c + 1)
                        if self.__lattice.get_Lattice().has_edge(i, j):
                            self.__lattice.get_Lattice().remove_edge(i, j)
                        if self.__lattice.get_Lattice().has_edge(j, i):
                            self.__lattice.get_Lattice().remove_edge(j, i)
                        self.__lattice.get_Lattice().node[new_pos]['pos'] = None
                        self.__lattice.get_Lattice().node[new_pos]['Contr'] += 1
                        self.__lattice.get_Lattice().node[j]['pos'] = None
            self.__lattice.get_Lattice().remove_node(i)
            self.__lattice.add_lostNode(i)

        return True

    # The function to call for executing a Fiber Walk simulation
    def walk(self, steps = 100, frequency = 2, position ='center'):

        # init all Walks

        for n in range(self.__objects):
            restart_count = 0
            step_possible = True
            # self.resetAllArrays()
            self.reset_arrays(True)
            self.__contr_count = 0
            self.__sa_path_count = 0
            self.__sa_count = 0
            seed_pos = self.define_seed(position)
            self.__walk.add_node(seed_pos, pos=(0., 0., 0.), counter=1)
            # arrays to collect statistics
            self.__positions.append([seed_pos, seed_pos, False])
            self.__all_distances_to_origin.append([0.0000001])
            self.__all_valid_neighbors.append([0])
            self.__all_contract_neighbors.append([0])
            self.__all_sa.append([0])
            self.__all_sa_paths.append([0])
            self.__all_contr.append([0])
            restart_active = False
            step_count_tmp = len(self.__distances_to_origin)
            step_count_max = 0
            old_step_count = 0
            # just one step downwards to init the walk
            # some init for the run
            new_walk = self.__positions[0][2]
            start_pos = self.__positions[0][0]
            self.__distances_to_origin = self.__all_distances_to_origin[n]
            self.__valid_neighbors = self.__all_valid_neighbors[n]
            self.__contract_neighbors = self.__all_contract_neighbors[n]
            self.__sa = self.__all_sa[n]
            self.__sa_path = self.__all_sa_paths[n]
            self.__contr = self.__all_contr[n]
            self.__sa_count = self.__sa[len(self.__sa) - 1]
            self.__sa_path_count = self.__sa_path[len(self.__sa_path) - 1]
            self.__contr_count = self.__contr[len(self.__contr) - 1]
            old_pos = start_pos

            while step_possible:

                # make array of currentPos to increment
                newPos = self.__positions[0][1]

                if not restart_active:
                    print 'step: ' + str(step_count_tmp) + ' of Walk Nr. ' + str(n + 1) + '/' + str(
                        self.__objects)
                    step_count_tmp += 1
                    if step_count_max - step_count_tmp > 20:
                        step_count_max = step_count_tmp
                restart_active = False
                # just the first step
                if newPos == seed_pos:
                    newPos, old_pos, stop, new_walk = self.init_seed_contraction(newPos, old_pos)
                    if len(self.__distances_to_origin) < steps and not stop:
                        step_possible = True
                    else:
                        step_possible = False
                    if stop:
                        step_possible = False
                        restart_active = True
                    self.__positions[0] = [old_pos, newPos, new_walk]
                # if step is not trapped, then old_pos and new_pos are different
                elif old_pos != newPos:
                    # save current state to recover it in case of a restart
                    nx.write_gpickle(self.__walk, './tmp/wtmp' + str(step_count_tmp))
                    nx.write_gpickle(self.__lattice, './tmp/ltmp' + str(step_count_tmp))
                    outfile = open('./tmp/WPAtmp' + str(step_count_tmp), 'wb')
                    pickle.dump(self.__positions[0], outfile)
                    outfile.close()
                    outfile = open('./tmp/nptmp' + str(step_count_tmp), 'wb')
                    pickle.dump(newPos, outfile)
                    outfile.close()
                    outfile = open('./tmp/optmp' + str(step_count_tmp), 'wb')
                    pickle.dump(old_pos, outfile)
                    outfile.close()
                    outfile = open('./tmp/nwtmp' + str(step_count_tmp), 'wb')
                    pickle.dump(new_walk, outfile)
                    outfile.close()
                    outfile = open('./tmp/dtotmp' + str(step_count_tmp), 'wb')
                    pickle.dump(self.__distances_to_origin, outfile)
                    outfile.close()
                    outfile = open('./tmp/satmp' + str(step_count_tmp), 'wb')
                    pickle.dump(self.__sa, outfile)
                    outfile.close()
                    outfile = open('./tmp/saptmp' + str(step_count_tmp), 'wb')
                    pickle.dump(self.__sa_path, outfile)
                    outfile.close()
                    outfile = open('./tmp/ctmp' + str(step_count_tmp), 'wb')
                    pickle.dump(self.__contr, outfile)
                    outfile.close()
                    newPos, old_pos, stop, new_walk = self.select_walk_edge_leading_label_self_avoid_contract_sq(newPos,
                                                                                                               old_pos,
                                                                                                               new_walk)
                elif old_pos == newPos:
                    print 'Walk reached stopping configuration at step ' + str(step_count_tmp)
                    # restart condition
                    if len(self.__distances_to_origin) < steps + 1:
                        if step_count_tmp > 2:
                            if restart_count < steps:
                                print 'restart the Fiber Walk'
                                if old_step_count >= step_count_tmp:
                                    step_count_tmp = step_count_max - 1
                                    # preventing getting negative step numbers. Step 1 is never trapped.
                                    if step_count_tmp == 1:
                                        step_count_tmp = 2
                                elif old_step_count < step_count_tmp:
                                    old_step_count = step_count_tmp
                                    # set the step counter one step back
                                    # Note: This occurs until the walk is longer then the first trapping
                                    step_count_tmp -= 1
                                    if step_count_tmp == 1:
                                        step_count_tmp = 2
                                # recover previous state
                                self.__walk = nx.read_gpickle('./tmp/wtmp' + str(step_count_tmp))
                                self.__lattice = nx.read_gpickle('./tmp/ltmp' + str(step_count_tmp))
                                infile = open('./tmp/WPAtmp' + str(step_count_tmp), 'rb')
                                self.__positions[0] = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/optmp' + str(step_count_tmp), 'rb')
                                old_pos = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/nptmp' + str(step_count_tmp), 'rb')
                                newPos = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/nwtmp' + str(step_count_tmp), 'rb')
                                new_walk = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/dtotmp' + str(step_count_tmp), 'rb')
                                self.__distances_to_origin = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/satmp' + str(step_count_tmp), 'rb')
                                self.__sa = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/saptmp' + str(step_count_tmp), 'rb')
                                self.__sa_path = nx.read_gpickle(infile)
                                infile.close()
                                infile = open('./tmp/ctmp' + str(step_count_tmp), 'rb')
                                self.__contr = nx.read_gpickle(infile)
                                infile.close()
                                restart_count = +1
                                restart_active = True
                                step_count_max = step_count_tmp
                    elif stop:
                        restart_active = True

                if len(self.__distances_to_origin) < steps + 1:
                    step_possible = True
                else:
                    step_possible = False

                if not restart_active:
                    # extend the current walk stats
                    self.__positions[0] = [old_pos, newPos, new_walk]
                    d = np.linalg.norm(np.array(newPos, dtype=float) - np.array(seed_pos, dtype=float), self.__norm)
                    d = d * d
                    self.__distances_to_origin.append(d)
                    self.__sa.append(self.__sa_count)
                    self.__sa_path.append(self.__sa_path_count)
                    self.__contr.append(self.__contr_count)
                    self.__positions[0] = [old_pos, newPos, new_walk]

                # store overall statistics
                self.__all_distances_to_origin[n] = self.__distances_to_origin
                self.__all_valid_neighbors[n] = self.__valid_neighbors
                self.__all_contract_neighbors[n] = self.__contract_neighbors
                self.__all_sa[n] = self.__sa
                self.__all_sa_paths[n] = self.__sa_path
                self.__all_contr[n] = self.__contr

            # extend stats arrays
            self.__all_walks_sa.append(self.__all_sa)
            self.__all_walks_sa_paths.append(self.__all_sa_paths)
            self.__allWalksContr.append(self.__all_contr)
            self.__all_walk_distances.append(self.__all_distances_to_origin)
