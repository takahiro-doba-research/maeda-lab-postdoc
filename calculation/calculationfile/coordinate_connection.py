import networkx as nx
import numpy as np

from .coordinate import Coordinate, concat


class CoordinateConnection(Coordinate):

    def __init__(self, df=None):
        super().__init__(df)

        """
        Examples
        --------
        self.my_name = 0
        self.other_names = [1, 2, 3, 4, 5]
        self.ref_labels = [
            (0, 1, {'ref_labels': [1, 2, 3]}),
            (0, 2, {'ref_labels': [20, 21, 22]}),
            (0, 3, {'ref_labels': [29, 31, 32]}),
            (0, 4, {'ref_labels': [36, 38, 39]}),
            (0, 5, {'ref_labels': [52, 47, 42]})
        ]
        self.drop_labels = [
            (0, 1, {'drop_labels': [2, 3, 4, 5, 6, 8, 9, 10, 11, 12]}),
            (0, 2, {'drop_labels': [21, 22, 23, 24]}),
            (0, 3, {'drop_labels': [31, 32, 33, 34]}),
            (0, 4, {'drop_labels': [38, 39, 40, 41]}),
            (0, 5, {'drop_labels': [42, 43, 44, 45, 46, 47, 48, 49, 50, 51]})
        ]
        self.bond_length = [
            (0, 1, {'bond_length': 1.4217396286293125}),
            (0, 2, {'bond_length': 1.4155479175515955}),
            (0, 3, {'bond_length': 1.5106539178425136}),
            (0, 4, {'bond_length': 1.4805055760779424}),
            (0, 5, {'bond_length': 1.3400573916193452})
        ]
        self.step = [
            (0, 1, {'step': 1}),
            (0, 2, {'step': 1}),
            (0, 3, {'step': 3}),
            (0, 4, {'step': 1}),
            (0, 5, {'step': 3})
        ]
        self.angles = [
            (0, 1, {'angles': [0]}),
            (0, 2, {'angles': [0]}),
            (0, 3, {'angles': [0, 120, 240]}),
            (0, 4, {'angles': [0]}),
            (0, 5, {'angles': [0, 120, 240]})
        ]
        """

    def read_connectivity(self):
        labels = list(self.df.index)
        notes = list(self.df["note"])
        lines = {
            note[0]:[note[1], label, *note[2:]]
            for label, note in zip(labels, notes) if len(note) >= 3
        }
        self.other_names = list(lines.keys())
        self.my_name = list(
            {note[0] for note in notes} - set(self.other_names)
        )[0]
        self.ref_labels = [
            (self.my_name, other_name, {"ref_labels":line[0:3]})
            for other_name, line in lines.items()
        ]
        self.drop_labels = [
            (self.my_name, other_name,
             {"drop_labels":[label for label, note in zip(labels, notes)
              if note[0] == other_name]})
            for other_name in self.other_names
        ]
        self.bond_length = [
            (self.my_name, other_name,
             {"bond_length":self.get_distance(line[0], line[1])})
            for other_name, line in lines.items()
        ]
        self.step = [
            (self.my_name, other_name, {"step":line[3]}) if len(line) == 4
            else (self.my_name, other_name, {"step":1})
            for other_name, line in lines.items()
        ]
        self.dihedral_angles = [
            (step[0], step[1],
             {"dihedral_angles":list(range(0, 360, int(360/step[2]["step"])))})
            for step in self.step
        ]


def connect(coordinates, distances=None, dihedral_angles=None):
    """
    Parameters
    ----------
    coordinates : list of CoordinateConnection object
        "read_connectivity" must be finished.
    distances : list of tuples
    dihedral_angles : list of tuples
    
    Examples
    --------
    distances = [
        (0, 1, {'distance': 1.4217396286293125}),
        (0, 2, {'distance': 1.4155479175515955}),
        (0, 3, {'distance': 1.5106539178425136}),
        (0, 4, {'distance': 1.4805055760779424}),
        (0, 5, {'distance': 1.3400573916193452})
    ]
    dihedral_angles = [
        (0, 1, {'dihedral_angle': 0}),
        (0, 2, {'dihedral_angle': 0}),
        (0, 3, {'dihedral_angle': 120}),
        (0, 4, {'dihedral_angle': 0}),
        (0, 5, {'dihedral_angle': 240})
    ]
    """
    # Make info0
    # edge attributes : "ref_labels", "drop_labels", "bond_length"
    info0 = nx.DiGraph()

    for coordinate in coordinates:
        info0.add_edges_from(coordinate.ref_labels)
        info0.add_edges_from(coordinate.drop_labels)
        info0.add_edges_from(coordinate.bond_length)

    # Make info1
    # edge attributes : "distance", "dihedral_angle"
    info1 = nx.Graph()
    info1.add_edges_from(info0.edges)

    for node0, node1 in info1.edges:
        node_s, node_l = sorted([node0, node1])
        info1.add_edge(
            node0, node1, distance=info0.edges[node_l, node_s]["bond_length"]
        )
        info1.add_edge(node0, node1, dihedral_angle=0)

    # Overwrite info1 with given distances and dihedral angles
    if distances != None:
        info1.add_edges_from(distances)
    if dihedral_angles != None:
        info1.add_edges_from(dihedral_angles)

    # Make a product graph
    graph = nx.Graph()
    graph.add_nodes_from(
        [(coordinate.my_name, {"coordinate":Coordinate(coordinate.df.copy())})
         for coordinate in coordinates]
    )

    # Connect the product graph using information of info0 and info1
    for node0, node1, attr in info1.edges(data=True):
        matrix0 = graph.nodes[node0]["coordinate"].get_to_x_axis(
            info0.edges[node0, node1]["ref_labels"],
            0,
            0
        )
        matrix1 = graph.nodes[node1]["coordinate"].get_to_x_axis_minus(
            info0.edges[node1, node0]["ref_labels"],
            attr["distance"],
            attr["dihedral_angle"]
        )
        graph.nodes[node0]["coordinate"].drop(
            info0.edges[node0, node1]["drop_labels"]
        )
        graph.nodes[node1]["coordinate"].drop(
            info0.edges[node1, node0]["drop_labels"]
        )
        
        component0 = nx.node_connected_component(graph, node0)
        component1 = nx.node_connected_component(graph, node1)
        
        for node in component0:
            graph.nodes[node]["coordinate"].transform(matrix0)
        for node in component1:
            graph.nodes[node]["coordinate"].transform(matrix1)
        
        graph.add_edge(node0, node1, **attr)

    coordinates = [attr["coordinate"] for _, attr in graph.nodes(data=True)]

    return concat(coordinates)

