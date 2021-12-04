from edge import QuadEdge


class Vertex():
    def __init__(self, coord, edge_ids:list, id):
        self.coord = coord
        self.edge_ids = edge_ids  # all 3-split created edges points towards new vertex
        self.id = id