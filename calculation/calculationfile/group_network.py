import copy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pyvis.network import Network
from scipy.stats import rankdata

from .eq_list import EQList
from .pt_list import PTList


class GroupNetwork:

    def __init__(self, name=None, rpn=None):
        self.name = name
        self.rpn = rpn
        self.graph = nx.Graph()

        for node0, nbrsdict in self.rpn.mdgraph.adjacency():
            eq0 = self.rpn.mdgraph.nodes[node0]["eq"]
            group0 = eq0.group

            if not self.graph.has_node(group0):
                self.graph.add_node(group0, eq_list=EQList(name=group0))
            
            self.graph.nodes[group0]["eq_list"][node0] = eq0

            for node1, keydict in nbrsdict.items():
                eq1 = self.rpn.mdgraph.nodes[node1]["eq"]
                group1 = eq1.group

                if not self.graph.has_node(group1):
                    self.graph.add_node(group1, eq_list=EQList(name=group1))

                self.graph.nodes[group1]["eq_list"][node1] = eq1

                for key, eattr in keydict.items():
                    pt = eattr["pt"]

                    if not self.graph.has_edge(group0, group1):
                        self.graph.add_edge(group0, group1, pt_list=PTList())

                    self.graph.edges[group0, group1]["pt_list"][pt.name] = pt

    def __getitem__(self, group):
        return self.graph.nodes[group]

    def remove_groups(self, groups):
        groups = copy.deepcopy(groups)
        self.graph.remove_nodes_from(groups)

    def shortest_path(self, source, target):
        return nx.shortest_path(self.graph, source, target)
        
    def set_traffic_max(self):

        for _, attr in self.graph.nodes(data=True):
            eq = max(attr["eq_list"].values(), key=lambda x:x.traffic)
            attr["traffic_max"] = eq.traffic

    def set_population_sum(self):

        for _, attr in self.graph.nodes(data=True):
            populations = [eq.population for eq in attr["eq_list"].values()]
            attr["population_sum"] = sum(populations)

    def set_energy(self):

        for _, nattr in self.graph.nodes(data=True):
            eq = min(nattr["eq_list"].values(), key=lambda x:x.E)
            nattr["energy"] = eq.E

        for _, _, eattr in self.graph.edges(data=True):
            pt = min(eattr["pt_list"].values(), key=lambda x:x.E)
            eattr["energy"] = pt.E

    def set_eq(self):

        for _, nattr in self.graph.nodes(data=True):
            eq = min(nattr["eq_list"].values(), key=lambda x:x.E)
            nattr["eq"] = eq

    def set_pt(self):

        for _, _, eattr in self.graph.edges(data=True):
            pt = min(eattr["pt_list"].values(), key=lambda x:x.E)
            eattr["pt"] = pt

    def adj_groups(self, groups):

        if isinstance(groups, int):
            groups = [groups]

        adj_groups = set()

        for node0, node1 in self.graph.edges(data=False):

            if node0 in groups or node1 in groups:
                adj_groups = adj_groups | {node0, node1}

        return sorted(list(adj_groups - set(groups)))

    def adj_groups_plus_one(self, groups):

        if isinstance(groups, int):
            groups = [groups]

        adj_groups = self.adj_groups(groups)
        adj_groups_plus_one = self.adj_groups(adj_groups)
        return sorted(list(set(adj_groups_plus_one) - set(groups)))

    def print_path(self, group0, group1, E=float("inf")):
        eqs0 = self.rpn.eqs_in_group(group0)
        eqs1 = self.rpn.eqs_in_group(group1)

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.rpn.mdgraph.edges(data="pt", keys=True):

            if (node0 in eqs0 and node1 in eqs1) or (node0 in eqs1 and node1 in eqs0) and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.rpn.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.rpn.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def show(self, path="default"):
        G = nx.Graph()
        n_nodes = self.graph.number_of_nodes()
        R = 70*n_nodes / (2*np.pi)

        for node, nattr in self.graph.nodes(data=True):
            cmap = plt.cm.jet(nattr["energy"]/200)
            cmap255 = (cmap[0]*255, cmap[1]*255, cmap[2]*255, 0.7)
            eq_names = sorted(nattr['eq_list'].keys())
            eq_names = "\n".join([
                ", ".join([
                    str(name) for name in eq_names[i*10:(i+1)*10]
                ])
                for i in range(len(eq_names)//10 + 1)
            ])
            T = 2*np.pi*(sorted(self.graph.nodes).index(node)/n_nodes - 1/4)
            G.add_node(
                node,
                label=str(node),
                shape="box",
                color=f"rgb{cmap255}",
                physics=False,
                title=f"EQ{eq_names}",
                x=R*np.cos(T),
                y=R*np.sin(T)
            )

        for node0, node1, eattr in self.graph.edges(data=True):
            pt_names = sorted(eattr['pt_list'].keys())
            pt_names = "\n".join([
                ", ".join([
                    str(name) for name in pt_names[i*10:(i+1)*10]
                ])
                for i in range(len(pt_names)//10 + 1)
            ])
            G.add_edge(
                node0,
                node1,
                physics=False,
                title=f"PT{pt_names}"
            )

        net = Network(height=800, width=1400)
        net.from_nx(G)

        if path == "default":
            path = f"{self.name}.html"

        net.show(path)

        # Remove the loading bar from the html file.
        with open(path, "r") as f:
            lines = f.readlines()

        index = [
            i for i, line in enumerate(lines)
            if '<div id="loadingBar">' in line
        ]

        if index:
            lines[index[0]:index[0]+9] = []

        with open(path, "w") as f:
            f.writelines(lines)

    def show_test(self, path):
        G = nx.Graph()
        n_nodes = self.graph.number_of_nodes()

        for node, nattr in self.graph.nodes(data=True):
            cmap = plt.cm.jet(nattr["energy"]/200)
            cmap255 = (cmap[0]*255, cmap[1]*255, cmap[2]*255, cmap[3])
            eq_names = sorted(nattr['eq_list'].keys())
            eq_names = "\n".join([
                ", ".join([
                    str(name) for name in eq_names[i*10:(i+1)*10]
                ])
                for i in range(len(eq_names)//10 + 1)
            ])
            G.add_node(
                node,
                label=str(node),
                shape="box",
                color=f"rgb{cmap255}",
                physics=True,
                title=f"EQ{eq_names}",
            )

        for node0, node1, eattr in self.graph.edges(data=True):
            pt_names = sorted(eattr['pt_list'].keys())
            pt_names = "\n".join([
                ", ".join([
                    str(name) for name in pt_names[i*10:(i+1)*10]
                ])
                for i in range(len(pt_names)//10 + 1)
            ])
            G.add_edge(
                node0,
                node1,
                physics=True,
                title=f"PT{pt_names}"
            )

        net = Network(height=800, width=1400)
        net.from_nx(G)
        net.show(path)

        """
        # Remove the loading bar from the html file.
        with open(path, "r") as f:
            lines = f.readlines()

        index = [
            i for i, line in enumerate(lines)
            if '<div id="loadingBar">' in line
        ]

        if index:
            lines[index[0]:index[0]+9] = []

        with open(path, "w") as f:
            f.writelines(lines)
        """