import copy
import itertools

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from networkx.exception import NetworkXError
from pyvis.network import Network
from scipy.stats import rankdata

from .eq_list import EQList
from .pt_list import PTList
from .eq_pref import EQPref
from .eq_popl import EQPopl
from .sim import Sim
from .eq_cl import EQCl
from .eq_pcl import EQPcl


class ReactionPathNetwork:

    def __init__(self, name=None, eq_list=None, pt_list=None,
            eq_pref=None, eq_popl=None, sim=None,
            eq_cl=None, eq_pcl=None):
        self.name = name
        self.eq_list = eq_list
        self.pt_list = pt_list

        self.mdgraph = nx.MultiDiGraph()
        self.mdgraph.add_nodes_from([
            (name, {"eq":eq}) 
            for name, eq in self.eq_list.items()
        ])
        self.mdgraph.add_edges_from([
            (pt.connection[0], pt.connection[1], {"pt":pt})
            for pt in self.pt_list.values()
        ])

        try:
            self.mdgraph.remove_node("??")
        except NetworkXError:
            # When there is no "??"
            pass

        self.eq_pref = eq_pref
        self.eq_popl = eq_popl
        self.sim = sim
        self.eq_cl = eq_cl
        self.eq_pcl = eq_pcl

        if isinstance(self.eq_pref, EQPref):
            mapping = self.eq_pref.to_dict()
            self.set_eq_pref(mapping)
        
        if isinstance(self.eq_popl, EQPopl):
            mapping = self.eq_popl.to_dict()
            self.set_eq_popl(mapping)

        if isinstance(self.sim, Sim):
            mapping = self.sim.to_dict()
            self.set_traffic(mapping)

        if isinstance(self.eq_cl, EQCl):
            mapping = self.eq_cl.to_dict()
            self.set_group(mapping)

        if isinstance(self.eq_pcl, EQPcl):
            mapping = self.eq_pcl.to_dict()
            self.set_group(mapping)

    def set_group(self, mapping):

        for name, group in mapping.items():
            self.mdgraph.nodes[name]["eq"].group = group

    def set_eq_pref(self, mapping):

        for name, preference in mapping.items():
            self.mdgraph.nodes[name]["eq"].preference = preference

        for _, eq in self.mdgraph.nodes(data="eq"):
            if not hasattr(eq, "preference"):
                eq.preference = 0.0

    def set_eq_popl(self, mapping):

        for name, population in mapping.items():
            self.mdgraph.nodes[name]["eq"].population = population

        for _, eq in self.mdgraph.nodes(data="eq"):
            if not hasattr(eq, "population"):
                eq.population = 0.0

    def set_traffic(self, mapping):

        for name, traffic in mapping.items():
            self.mdgraph.nodes[name]["eq"].traffic = traffic

        for _, eq in self.mdgraph.nodes(data="eq"):
            if not hasattr(eq, "traffic"):
                eq.traffic = 0.0

    def adjust_energy(self, name=None):
        """
        Adjusts energies to kJ/mol.
        """
        if name == None:
            E_std = min(
                self.mdgraph.nodes(data="eq"), key=lambda x:x[1].E
            )[1].E

        else:
            E_std = self.mdgraph.nodes[name]["eq"].E

        for node in self.mdgraph.nodes:
            self.mdgraph.nodes[node]["eq"].E = (
                self.mdgraph.nodes[node]["eq"].E - E_std
            ) * 2625.5

        for edge in self.mdgraph.edges:
            self.mdgraph.edges[edge]["pt"].E = (
                self.mdgraph.edges[edge]["pt"].E - E_std
            ) * 2625.5

    def most_stable_eq(self):
        return min(self.mdgraph.nodes(data="eq"), key=lambda x:x[1].E)[0]

    def energy_of_eq(self, name):
        return self.mdgraph.nodes[name]["eq"].E

    def group_of_eq(self, name):
        return self.mdgraph.nodes[name]["eq"].group

    def preference_of_eq(self, name):
        return self.mdgraph.nodes[name]["eq"].preference

    def population_of_eq(self, name):
        return self.mdgraph.nodes[name]["eq"].population

    def traffic_of_eq(self, name):
        return self.mdgraph.nodes[name]["eq"].traffic

    def eqs_in_group(self, group):
        names = [
            node for node, eq in self.mdgraph.nodes(data="eq")
            if eq.group == group
        ]
        return names

    def adj_eqs(self, names, E=float("inf")):

        if isinstance(names, int):
            names = [names]

        nodes = []

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node0 in names or node1 in names and pt.E < E:
                nodes += [node0, node1]

        nodes = sorted(list(set(nodes) - set(names)))

        return nodes

    def pred_eqs(self, names, E=float("inf")):

        if isinstance(names, int):
            names = [names]

        nodes = []

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node1 in names and pt.E < E:
                nodes += [node0, node1]

        nodes = sorted(list(set(nodes) - set(names)))

        return nodes

    def succ_eqs(self, names, E=float("inf")):

        if isinstance(names, int):
            names = [names]

        nodes = []

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node0 in names and pt.E < E:
                nodes += [node0, node1]

        nodes = sorted(list(set(nodes) - set(names)))

        return nodes
    
    def adj_groups(self, groups, E=float("inf")):

        if isinstance(groups, int):
            groups = [groups]

        adj_groups = []

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if (self.group_of_eq(node0) in groups
                    or self.group_of_eq(node1) in groups
                    and pt.E < E):
                adj_groups += [
                    self.group_of_eq(node0), self.group_of_eq(node1)
                ]

        adj_groups = sorted(list(set(adj_groups) - set(groups)))

        return adj_groups

    def shortest_path(self, source, target):
        mgraph = self.mdgraph.to_undirected()
        return nx.shortest_path(mgraph, source, target)

    def shortest_path_plus_one(self, source, target):
        eqs = self.shortest_path(source, target)
        eqs_plus_one = []

        for node0, node1, _ in self.mdgraph.edges(data=False, keys=True):
            if node0 in eqs or node1 in eqs:
                eqs_plus_one += [node0, node1]
        
        return sorted(list(set(eqs_plus_one)))

    def print_summary(self):
        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):
            path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                node0, pt.name, node1
            )
            energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                self.mdgraph.nodes[node0]["eq"].E,
                pt.E,
                self.mdgraph.nodes[node1]["eq"].E
            )
            print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def print_adj(self, names, E=float("inf")):
        
        if isinstance(names, int):
            names = [names]

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node0 in names or node1 in names and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def print_pred(self, names, E=float("inf")):

        if isinstance(names, int):
            names = [names]

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node1 in names and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def print_succ(self, names, E=float("inf")):
        
        if isinstance(names, int):
            names = [names]

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node0 in names and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def print_path(self, eqs0, eqs1, E=float("inf")):

        if isinstance(eqs0, int):
            eqs0 = [eqs0]
        if isinstance(eqs1, int):
            eqs1 = [eqs1]

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if node0 in eqs0 and node1 in eqs1 and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")

    def print_pt_path(self, pts, E=float("inf")):

        if isinstance(pts, int):
            pts = [pts]

        print("  EQ ->   PT ->   EQ  :   E_EQ ->  E_PT ->  E_EQ")
        print("------------------------------------------------")

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):

            if pt.name in pts and pt.E < E:
                path = "{: >4d} -> {: >4d} -> {: >4d}".format(
                    node0, pt.name, node1
                )
                energy = "{: >5.1f} -> {: >5.1f} -> {: >5.1f}".format(
                    self.mdgraph.nodes[node0]["eq"].E,
                    pt.E,
                    self.mdgraph.nodes[node1]["eq"].E
                )
                print(f"{path:20s}  :  {energy}")

        print("------------------------------------------------")
        print("                   (electronic energy in kJ/mol)")
    
    def to_gml(self, path):
        output_graph = copy.deepcopy(self.mdgraph)

        for _, eq in output_graph.nodes(data="eq"):
            del eq["coordinate"], eq["normal_mode"]

        for _, _, _, pt in output_graph.edges(data="pt", keys=True):
            del pt["coordinate"], pt["normal_mode"]

        nx.write_gml(output_graph, path)

    def show(self, path="default"):
        G = nx.MultiDiGraph()
        n_nodes = self.mdgraph.number_of_nodes()
        R = 70*n_nodes / (2*np.pi)

        for node, eq in self.mdgraph.nodes(data="eq"):
            cmap = plt.cm.jet(eq.E/200)
            cmap255 = (cmap[0]*255, cmap[1]*255, cmap[2]*255, cmap[3])
            T = 2*np.pi*(node/n_nodes - 1/4)

            G.add_node(
                node,
                label=str(node),
                shape="circle",
                color=f"rgb{cmap255}",
                physics=False,
                title=f"EQ{eq.name}",
                x=R*np.cos(T),
                y=R*np.sin(T)
            )

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):
            G.add_edge(
                node0,
                node1,
                physics=False,
                title=f"PT{pt.name}"
            )

        net = Network(height=800, width=1400, directed=True)
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

    def show2(self, path="default"):
        G = nx.MultiDiGraph()
        #group_eqs = self.eq_pcl.to_group_dict()
        group_eqs = {}

        for _, eq in self.mdgraph.nodes(data="eq"):
            group = eq.group
            
            if group in group_eqs.keys():
                group_eqs[group].append(eq.name)

            else:
                group_eqs[group] = [eq.name]

        group_num = {group:len(eqs) for group, eqs in group_eqs.items()}
        group_r = {group:70*num/(2*np.pi) for group, num in group_num.items()}
        group_r_cumsum = {
            group:sum([v for k, v in group_r.items() if k <= group])
            for group in group_r.keys()
        }
        R = 3.0 * max(group_r_cumsum.values()) / np.pi

        for node, eq in self.mdgraph.nodes(data="eq"):
            cmap = plt.cm.jet(eq.E/200)
            cmap255 = (cmap[0]*255, cmap[1]*255, cmap[2]*255, cmap[3])
            group = eq.group

            if group_num[group] == 1:
                r = 0
            else:
                r = group_r[group]

            T = 2*np.pi*(
                3.0*(2*group_r_cumsum[group] - group_r[group])
                /(2*np.pi*R) - 1/4
            )
            t = 2*np.pi*(
                sorted(group_eqs[group]).index(eq.name)
                /group_num[group] - 1/4
            )
            
            G.add_node(
                node,
                label=str(node),
                shape="dot",
                color=f"rgb{cmap255}",
                physics=False,
                title=f"G{group}",
                x=R*np.cos(T)+r*np.cos(t),
                y=R*np.sin(T)+r*np.sin(t)
            )

        for node0, node1, _, pt in self.mdgraph.edges(data="pt", keys=True):
            G.add_edge(
                node0,
                node1,
                physics=False,
                title=f"PT{pt.name}"
            )

        net = Network(height=800, width=1400, directed=True)
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
