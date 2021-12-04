from edge import QuadEdge


class Face():
    def __init__(self, edge_id_rot:list, id):
        self.edge_id_rot = edge_id_rot  # tuple (edge id, rot num), ordered CCW
        self.id = id
