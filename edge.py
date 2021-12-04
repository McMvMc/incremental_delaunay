SYM = 2
N_ROT = 4


def check_rot(left_id, src_e):
    rot = [src_e.left_id, src_e.orig_id, src_e.right_id, src_e.dest_id].index(left_id)
    return rot


class QuadEdge:
    def __init__(self, orig_id, dest_id, left_id, right_id, id):
        self.orig_id = orig_id
        self.dest_id = dest_id
        self.left_id = left_id
        self.right_id = right_id

        self.id = id

    def get_endpt_ids(self, rot_num):
        tmp = [[self.orig_id, self.dest_id],
               [self.right_id, self.left_id],
               [self.dest_id, self.orig_id],
               [self.left_id, self.right_id]]
        return tmp[rot_num]

    def get_right_id(self, rot_num):
        tmp = [self.right_id, self.dest_id, self.left_id, self.orig_id]
        return tmp[rot_num]

    def get_left_id(self, rot_num):
        tmp = [self.left_id, self.orig_id, self.right_id, self.dest_id]
        return tmp[rot_num]

    def set_left(self, id, rot_num):
        if rot_num == 0:
            self.left_id = id
        elif rot_num == 1:
            self.orig_id = id
        elif rot_num == 2:
            self.right_id = id
        elif rot_num == 3:
            self.dest_id = id

    def set_right(self, id, rot_num):
        if rot_num == 0:
            self.right_id = id
        elif rot_num == 1:
            self.dest_id = id
        elif rot_num == 2:
            self.left_id = id
        elif rot_num == 3:
            self.orig_id = id

    # Note: CCW operators
    def Rnext(self, all_tri, rot_num):
        "next edge around right face, with same right face"
        cur_right_id = self.get_right_id(rot_num)
        cur_idx = all_tri[cur_right_id].edge_ids.index(self.id)
        rnext_idx = all_tri[cur_right_id].edge_ids[(cur_idx+1)%3]
        return rnext_idx

    def Lnext(self, all_tri, all_edges, rot_num):
        "next edge around left face, with same left face"
        cur_left_id = self.get_left_id(rot_num)
        try:
            for i, e_id_rot in enumerate(all_tri[cur_left_id].edge_id_rot):
                if e_id_rot[0] == self.id:
                    cur_idx = i
            lnext_id = all_tri[cur_left_id].edge_id_rot[(cur_idx+1)%3][0]
            out_e = all_edges[lnext_id]
            try:
                out_rot = check_rot(cur_left_id, out_e)
                return out_e, out_rot
            except:
                raise Exception(f"failed to determine rot for edge {out_e.id}")
        except:
            raise Exception(f"Lnext not found for edge {self.id}")

    def Onext(self, all_tri, all_edges, rot_num):
        "next edge around origin, with same origin"
        cur_left_id = self.get_left_id(rot_num)
        cur_idx = all_tri[cur_left_id].edge_ids.index(self.id)
        lnext_idx = all_tri[cur_left_id].edge_ids[(cur_idx-1)]
        return lnext_idx

    def Dnext(self, all_tri, rot_num):
        "next edge around dest, with same dest"
        cur_right_id = self.get_right_id(rot_num)
        cur_idx = all_tri[cur_right_id].edge_ids.index(self.id)
        rnext_idx = all_tri[cur_right_id].edge_ids[(cur_idx-1)]
        return rnext_idx

    # Note: CW operators

    def Rprev(self, all_tri, rot_num):
        "prev edge around right face, with same right face"
        return self.Dnext(all_tri, rot_num)

    def Lprev(self, all_tri, rot_num):
        "prev edge around left face, with same left face"
        return self.Onext(all_tri, rot_num)

    def Oprev(self, all_tri, rot_num):
        "prev edge around origin, with same origin"
        return self.Rnext(all_tri, rot_num)

    def Dprev(self, all_tri, rot_num):
        "prev edge around dest, with same dest"
        return self.Lnext(all_tri, rot_num)
