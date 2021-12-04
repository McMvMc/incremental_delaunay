import numpy as np

import matplotlib.pyplot as plt

from edge import QuadEdge, SYM, N_ROT, check_rot
from face import Face
from vertex import Vertex

MAX = 30000.
MIN = -30000.


def check_left(orig:Vertex, dest:Vertex, pt:Vertex):
    det = np.array([[1., orig.coord[0], orig.coord[1]],
                    [1., dest.coord[0], dest.coord[1]],
                    [1., pt.coord[0], pt.coord[1]]])
    return np.linalg.det(det) > 0.


def check_in_tri(e0:list, e1:list, e2:list, pt:Vertex):
    return check_left(e0[0], e0[1], pt) and check_left(e1[0], e1[1], pt) and check_left(e2[0], e2[1], pt)


def draw_vert(vert, c):
    ax = plt.gca()
    vert_circ = plt.Circle((vert.coord[0], vert.coord[1]), 0.5, color=c, fill=True)
    ax.add_patch(vert_circ)
    plt.ion()
    plt.xlim([DISPLAY_MIN, DISPLAY_MAX])
    plt.ylim([DISPLAY_MIN, DISPLAY_MAX])
    plt.show()
    return vert_circ


def draw_circle(ccenter1, cradius1):
    ax = plt.gca()
    circle_near = plt.Circle((ccenter1[0], ccenter1[1]), cradius1 ** 0.5, color='b', fill=False)
    ax.add_patch(circle_near)
    plt.ion()
    plt.xlim([DISPLAY_MIN, DISPLAY_MAX])
    plt.ylim([DISPLAY_MIN, DISPLAY_MAX])
    plt.show()
    return circle_near


class Delaunay:
    def __init__(self):
        self.tri_list = {0: Face(edge_id_rot=[(0, 0), (1, 0), (2, 0)], id=0)}
        self.removed_tri_list = {}
        self.vert_list = {0: Vertex([MAX, MAX], edge_ids=[0,2], id=0),
                     1: Vertex([-MAX, MAX], edge_ids=[0,1], id=1),
                     2: Vertex([0, -MAX], edge_ids=[1,2], id=2)}
        self.edge_list = {0: QuadEdge(orig_id=0, dest_id=1, left_id=0, right_id=None, id=0),
                          1: QuadEdge(orig_id=1, dest_id=2, left_id=0, right_id=None, id=1),
                          2: QuadEdge(orig_id=2, dest_id=0, left_id=0, right_id=None, id=2)}

        self.max_face_id = 0
        self.max_edge_id = 2
        self.max_vert_id = 2

    def draw(self, draw_from_face=True):
        # l = []
        plt.clf()
        for tri_id in self.tri_list:
            edge_id_rot = self.tri_list[tri_id].edge_id_rot
            for e_id, e_rot in edge_id_rot:
                orig_id, dest_id = self.edge_list[e_id].get_endpt_ids(e_rot)
                if orig_id in [0,1,2] or dest_id in [0,1,2]:
                    continue
                x = self.vert_list[orig_id].coord[0]
                y = self.vert_list[orig_id].coord[1]
                dx = self.vert_list[dest_id].coord[0] - x
                dy = self.vert_list[dest_id].coord[1] - y
                plt.scatter(x, y, s=1.)
                plt.arrow(x=x, y=y, dx=dx, dy=dy, facecolor='red', width=0.5)
                # plt.annotate(f"v:{orig_id}", xy=(x,y))
                # plt.annotate(f"f:{tri_id} e:{e_id}", xy=(x+dx/1.3, y+dy/1.3))
        plt.ion()
        plt.xlim([DISPLAY_MIN, DISPLAY_MAX])
        plt.ylim([DISPLAY_MIN, DISPLAY_MAX])
        plt.show()

    def locate_pt(self, x):
        for face_idx in self.tri_list:
            e0_id_rot = self.tri_list[face_idx].edge_id_rot[0]
            e1_id_rot = self.tri_list[face_idx].edge_id_rot[1]
            e2_id_rot = self.tri_list[face_idx].edge_id_rot[2]

            orig_id0, dest_id0 = self.edge_list[e0_id_rot[0]].get_endpt_ids(e0_id_rot[1])
            orig_id1, dest_id1 = self.edge_list[e1_id_rot[0]].get_endpt_ids(e1_id_rot[1])
            orig_id2, dest_id2 = self.edge_list[e2_id_rot[0]].get_endpt_ids(e2_id_rot[1])
            o0, d0 = self.vert_list[orig_id0], self.vert_list[dest_id0]
            o1, d1 = self.vert_list[orig_id1], self.vert_list[dest_id1]
            o2, d2 = self.vert_list[orig_id2], self.vert_list[dest_id2]
            if check_in_tri([o0, d0], [o1, d1], [o2, d2], x):
                return face_idx
        return -1

    def remove_face(self, removed_tri):
        self.removed_tri_list[removed_tri.id] = removed_tri
        self.tri_list.pop(removed_tri.id)

    def split_to_3_tri(self, old_tri_id, x:Vertex):
        "create 3 new faces and 3 new edges. edges point towards x. edge id increase CCW"

        old_tri = self.tri_list[old_tri_id]

        # make 3 new faces and 3 new edges pointing towards x
        new_face_ids = [self.max_face_id + 1, self.max_face_id + 2, self.max_face_id + 3]
        new_edge_ids = [self.max_edge_id + 1, self.max_edge_id + 2, self.max_edge_id + 3]
        for new_idx, cur_edge_id_rot in enumerate(old_tri.edge_id_rot):
            cur_qedge = self.edge_list[cur_edge_id_rot[0]]
            cur_qedge_rot = 0 if cur_qedge.left_id == old_tri.id else SYM

            # new face
            new_e_id_rot_0 = (cur_qedge.id, cur_qedge_rot)
            new_e_id_rot_1 = (new_edge_ids[new_idx], 0)
            new_e_id_rot_2 = (new_edge_ids[new_idx-1], SYM)  # new edges pt towards x

            self.tri_list[new_face_ids[new_idx]] = Face(edge_id_rot=[new_e_id_rot_0, new_e_id_rot_1, new_e_id_rot_2],
                                                        id=new_face_ids[new_idx])
            cur_qedge.set_left(new_face_ids[new_idx], cur_qedge_rot)
            self.max_face_id = new_face_ids[new_idx]

            # new edge
            old_orig_id, old_dest_id = cur_qedge.get_endpt_ids(cur_edge_id_rot[1])
            self.edge_list[new_edge_ids[new_idx]] = QuadEdge(orig_id=old_dest_id, dest_id=x.id,
                                                             left_id=new_face_ids[new_idx],
                                                             right_id=new_face_ids[(new_idx+1)%3],
                                                             id=new_edge_ids[new_idx])
            self.max_edge_id += 1

            # add the new edge to the 3 vertices' edges list
            self.vert_list[old_dest_id].edge_ids.append(new_edge_ids[new_idx])
        x.edge_ids = new_edge_ids

        # remove from current tri list
        self.remove_face(old_tri)

        # self.draw()
        return

    def empty_circle(self, tri: Face, diag_vert:Vertex):
        # circum center
        e_id_rot0 = tri.edge_id_rot[0]
        e_id_rot1 = tri.edge_id_rot[1]

        e0_endpt_ids = self.edge_list[e_id_rot0[0]].get_endpt_ids(e_id_rot0[1])
        A = np.array(self.vert_list[e0_endpt_ids[0]].coord)
        B = np.array(self.vert_list[e0_endpt_ids[1]].coord)

        e1_endpt_ids = self.edge_list[e_id_rot1[0]].get_endpt_ids(e_id_rot1[1])
        assert e1_endpt_ids[0] == e0_endpt_ids[1], "end pt not in order for circumcenter calculation"
        C =np.array(self.vert_list[e1_endpt_ids[1]].coord)

        a2 = (B[0]-C[0]) ** 2 + (B[1]-C[1]) ** 2
        b2 = (A[0]-C[0]) ** 2 + (A[1]-C[1]) ** 2
        c2 = (A[0]-B[0]) ** 2 + (A[1]-B[1]) ** 2
        wa = a2 * (b2+c2-a2)
        wb = b2 * (a2+c2-b2)
        wc = c2 * (a2+b2-c2)

        ccenter = (wa * A + wb * B + wc * C) / (wa + wb + wc)

        # radius
        cradius2 = (A - ccenter) @ (A - ccenter).T

        dist = (np.array(diag_vert.coord) - ccenter) @ (np.array(diag_vert.coord) - ccenter).T

        is_in = np.linalg.det(np.array([[A[0], A[1], A[0]**2+A[1]**2, 1],
                                        [B[0], B[1], B[0]**2+B[1]**2, 1],
                                        [C[0], C[1], C[0]**2+C[1]**2, 1],
                                        [diag_vert.coord[0], diag_vert.coord[1],
                                         diag_vert.coord[0]**2+diag_vert.coord[1]**2, 1]]))

        # return dist >= cradius2, ccenter, cradius2
        return not is_in>=0, ccenter, cradius2

    def flip(self, tri_near, tri_far, diag_e, diag_e_rot, x, diag_vert, start_edge, far_eup, far_eup_rot):
        new_face_ids = [self.max_face_id+1, self.max_face_id+2]

        # all edge info
        start_edge_rot = check_rot(tri_near.id, start_edge)
        far_elo, far_elo_rot = far_eup.Lnext(self.tri_list, self.edge_list, rot_num=far_eup_rot)
        near_elo, near_elo_rot = diag_e.Lnext(self.tri_list, self.edge_list, rot_num=diag_e_rot)

        # add new diag
        self.max_edge_id += 1
        new_diag_edge = QuadEdge(orig_id=diag_vert.id, dest_id=x.id,
                                 left_id=new_face_ids[0], right_id=new_face_ids[1], id=self.max_edge_id)
        self.edge_list[new_diag_edge.id] = new_diag_edge

        # add upper diag face
        new_face_up = Face(edge_id_rot=[(start_edge.id, start_edge_rot),
                                        (far_eup.id, far_eup_rot),
                                        (new_diag_edge.id, 0)], id=new_face_ids[0])
        self.max_face_id = new_face_ids[0]
        self.tri_list[new_face_up.id] = new_face_up

        # add lower diag face
        new_face_lo = Face(edge_id_rot=[(new_diag_edge.id, 2),
                                        (far_elo.id, far_elo_rot),
                                        (near_elo.id, near_elo_rot)], id=new_face_ids[1])
        self.max_face_id = new_face_ids[1]
        self.tri_list[new_face_lo.id] = new_face_lo

        # remove 2 old faces
        self.remove_face(tri_near)
        self.remove_face(tri_far)

        # remove old diag edge from edge list and its verts' edge lists
        # try:
        self.vert_list[diag_e.orig_id].edge_ids.remove(diag_e.id)
        self.vert_list[diag_e.dest_id].edge_ids.remove(diag_e.id)
        # except:
        #     pass
        self.edge_list.pop(diag_e.id)

        # add new edge to x's edge list
        x.edge_ids.append(new_diag_edge.id)
        diag_vert.edge_ids.append(new_diag_edge.id)

        # change left of quadralateral edges to new faces
        start_edge.set_left(new_face_up.id, start_edge_rot)
        far_eup.set_left(new_face_up.id, far_eup_rot)
        far_elo.set_left(new_face_lo.id, far_elo_rot)
        near_elo.set_left(new_face_lo.id, near_elo_rot)

        return

    def find_not_delaunay(self, x: Vertex):
        was_delaunay = True
        # draw_vert(x, c='b')
        for e_id in x.edge_ids:
            cur_edge = self.edge_list[e_id]

            # get near far triangles
            tri_near_id = cur_edge.get_left_id(rot_num=SYM)
            tri_near = self.tri_list[tri_near_id]
            diag_e, diag_e_rot = cur_edge.Lnext(self.tri_list, self.edge_list, rot_num=SYM)
            tri_far_id = diag_e.get_right_id(diag_e_rot)
            if tri_far_id is None:
                continue
            tri_far = self.tri_list[tri_far_id]

            # upper left edge of the quadrilateral
            far_e0, far_e0_rot = diag_e.Lnext(self.tri_list, self.edge_list, rot_num=(SYM + diag_e_rot) % N_ROT)

            # lower left edge of the quadrilateral
            # far_e1, far_e1_rot = diag_e.Lnext(self.tri_list, self.edge_list, rot_num=(SYM + far_e0_rot) % N_ROT)

            # get vertices of the quadrilateral in CCW starting from new vert x
            diag_vert_id = far_e0.get_endpt_ids(rot_num=far_e0_rot)[1]
            diag_vert = self.vert_list[diag_vert_id]

            # test emty circle
            is_near_empty, ccenter1, cradius1 = self.empty_circle(tri_near, diag_vert)
            is_far_empty, ccenter2, cradius2 = self.empty_circle(tri_far, x)

            # circle_patch1 = draw_circle(ccenter1, cradius1)
            # circle_patch2 = draw_circle(ccenter2, cradius2)
            # circle_patch1.set_visible(False)
            # circle_patch2.set_visible(False)
            if not is_near_empty or not is_far_empty:
                self.flip(tri_near, tri_far, diag_e, diag_e_rot, x, diag_vert, cur_edge, far_e0, far_e0_rot)
                return False

        return was_delaunay

    def add_pt(self, x_coord):
        self.max_vert_id += 1
        x = Vertex(x_coord, edge_ids=[], id=self.max_vert_id)
        self.vert_list[self.max_vert_id] = x

        cur_tri = self.locate_pt(x)
        self.split_to_3_tri(cur_tri, x)
        # d.draw()
        # draw_vert(x, c='b')
        while True:
            was_delaunay = self.find_not_delaunay(x)
            if was_delaunay:
                break

        # draw_vert(x, c='g')

    def print_vert_edges(self):
        for v in self.vert_list:
            print(f"vert {v}: {self.vert_list[v].edge_ids}")

    def print_faces(self):
        for tri_id in self.tri_list:
            print(f"vert {tri_id}: {self.tri_list[tri_id].id}")


# MAX = 200.
# MIN = -200.
DISPLAY_MAX = 100
DISPLAY_MIN = -100
click_x, click_y, clicked = None, None, False
d = Delaunay()

# d.add_pt([4, -16])
# d.draw()
# d.add_pt([44, 17])
# d.draw()
# d.add_pt([-33, -24])
# d.draw()
# d.add_pt([-42, -23])
# d.draw()
# d.add_pt([30, 23])
# d.draw()


# MAX = 30.
# MIN = -30.
def onclick(event):
    # d.draw()
    global click_x, click_y, clicked
    click_x, click_y = event.xdata, event.ydata
    print('x = %d, y = %d' % (click_x, click_y))

    d.add_pt([click_x, click_y])
    d.draw()


fig, ax = plt.subplots()
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.xlim([DISPLAY_MIN, DISPLAY_MAX])
plt.ylim([DISPLAY_MIN, DISPLAY_MAX])
plt.show()
