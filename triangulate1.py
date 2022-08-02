"""
    File:        triangulate.py
    Author:      Madison Monroe
    Course:      CS 307 - Computational Geometry
    Assignment:  Problem Set 4 - Polygon Triangulation
    Description: Methods to generate y-monotone polygons,
    triangulate them, and verify the triangulation.
"""
from typing import List, Tuple

def order_vertices(polygon):
    """Given a list of vertices in clockwise order, creates a list of the
       indexes of polygon so that the vertices of polygon are in decreasing y
       order.

       Ex. polygon = [(-1,-3), (-2,-2), (0,0), (2,-4)]
       ordered = [2, 1, 0, 3]
       where [polygon[2], polygon[1], polygon[0], polygon[3]] is the list of
       coordinate ordered in decreasing-y order
       """
    top_index = bottom_index = left = right = None
    ordered_indexes = []

    if len(polygon) < 3:
        print("Error")
        return False

    index = 0
    #finds the vertex with the highest y coordinate
    while True:

        if polygon[index][1] < polygon[(index + 1) % len(polygon)][1]:
            index = (index + 1) % len(polygon)

        elif polygon[index][1] < polygon[(index - 1) % len(polygon)][1]:
            index = polygon[(index - 1) % len(polygon)]

        else:
            #top_index is index of the top vertex in polygon
            top_index = index
            print("top")
            left_index = (index - 1) % len(polygon)
            right_index = (index + 1) % len(polygon)
            ordered_indexes.append(index)
            break

    """traversing the left and right chains and inserting the higher vertex
       index into the order vertice list"""
    while polygon[left_index] != polygon[right_index]:
        if polygon[left_index][1] > polygon[right_index][1]:
            ordered_indexes.append(left_index)
            left_index = (left_index - 1) % len(polygon)

        else:
            ordered_indexes.append(right_index)
            right_index = (right_index + 1) % len(polygon)

    #adding the lowest vertice to the list
    ordered_indexes.append(left_index)
    #top_index is index of the lowest vertex in polygon
    bottom_index = left_index

    return (ordered_indexes, top_index, bottom_index)

def reflex(polygon, curr_index, stack, top_index, bottom_index):
    if len(stack) >= 2:

        curr_side = chain_side(curr_index, top_index, bottom_index)
        if chain_side(stack[-2], top_index, bottom_index) != curr_side \
           and chain_side(stack[-2], top_index, bottom_index) != 'both':
            return False

        elif orient(polygon[stack[-2]], polygon[stack[-1]], polygon[curr_index]) \
             == -1 and curr_side == 'left':
             return True

        elif orient(polygon[stack[-2]], polygon[stack[-1]], polygon[curr_index]) \
             == 1 and curr_side == 'right':
             return True

    if len(stack) == 1:
        return True

    return False

def chain_side(index, top_index, bottom_index):
    """Given the an index, and the indices of of the coordinate with the highest
       and lowest y-values, determine what side of the chain the index is on.
    """

    if top_index == index:
        return 'both'

    if top_index > bottom_index:
        if index > bottom_index and index < top_index:
            return 'left'

        return 'right'

    else: #bottom_index > top index
        if index > top_index and index < bottom_index:
            return 'right'

        return 'left'

def leftmost(polygon):
    """given a list of coordinates, find the minimum x_coordinate and the index of \
    this point"""
    minx = polygon[0][0]
    miny = polygon[0][1]
    minindex = 0

    for index in range(len(polygon)):
        x_coordinate, y_coordinate = polygon[index]
        if x_coordinate < minx:
            minx = x_coordinate
            miny = y_coordinate
            minindex = index

    return (((minx, miny), minindex))

def orient(p, q, r):
    """Determines whether point r is on the left, right, or collinear with p and
       q. Returns 1 if the point is to the left, -1 if it is to the right,
       and 0 if the point is collinear"""
    #splits tuple of points to get x and y coordinates
    px, py = p
    qx, qy = q
    rx, ry = r

    #finds the determinate as was done in class
    determinate = qx * ry + px * qy + rx * py - qx * py - rx * qy - ry * px

    if determinate > 0:
        return 1

    elif determinate < 0:
        return -1

    return 0

def find_slope(p, q):
    """Determines the slope given two points"""

    px, py = p
    qx, qy = q

    if px - qx == 0:
        return None

    return (py - qy) / (px - qx)

def intersection_found(linesegment1, linesegment2):
    """Determines if there is an intersections between two line segments"""
    point1 = linesegment1[0]
    point2 = linesegment1[1]
    point3 = linesegment2[0]
    point4 = linesegment2[1]

    #print(point1)
    #print(point2)
    #print(point3)
    #print(point4)
    #print(orient(point1, point2, point3))
    #print(orient(point1, point2, point4))
    #print(orient(point3, point4, point1))
    #print(orient(point3, point4, point2))
    orient1 = orient(point1, point2, point3)
    orient2 = orient(point1, point2, point4)
    orient3 = orient(point3, point4, point1)
    orient4 = orient(point3, point4, point2)

    if point1 == point3 and point2 == point4:
       return True

    if orient1 != orient2 and orient1 != 0 and orient2 != 0 \
      and (orient3 != orient4 and orient3 != 0 and orient4 != 0):
        return True

    return False

def findz(polygon, u, w, v):
    "Finds the z given a u, w, and v values. If there is no z, None is returned"

    intersection = False
    for index in range(len(polygon) - 1):
        if intersection_found((polygon[index], polygon[index + 1]), (u, w)):
            intersection = True
            break

        if not intersection:
            #print("no intersection")
            return (None, None)

    slope = find_slope(u, w)
    leftmost = None
    leftmost_index = 0

    for index in range(len(polygon)):
        #three calls to orient to determine if the vertex is in triange uvw
        vertex = polygon[index]

        if orient(v, w, vertex) == -1 and orient(w, u, vertex) == -1 and \
           orient(u, v, vertex) == -1:
            "finding another point on the line with the current vertice and \
            the same slope as line segment uw"
            x_coordinate, y_coordinate = vertex
            # b from y = mx + b
            b = y_coordinate - (slope * x_coordinate)
            #determining a new point that is above the vertex in the triangle
            #the line created by the vertex and the new point has the same slope
            #as wu so that it is as if we are sweeping a line parallel to uw
            newpointy = y_coordinate + 1
            newpointx = (newpointy - b) / slope
            "if the min is to the right of the line segment defined by \
            the vertex & (newpointx, newpointy) then we need to update the min"
            #print("here")
            #print(vertex)
            #print((newpointx, newpointy))
            #print("leftmost")
            #print(leftmost)
            if leftmost is None or orient(vertex, (newpointx, newpointy), leftmost) == -1:
                leftmost = vertex
                leftmost_index = index

    return (leftmost, leftmost_index)

def generate_monotone_polygon(num_vertices: int) -> List[Tuple[int,int]]:
    """
    Generate a y-monotone polygon on num_vertices vertices, returned
    as a list of points in clockwise order around its boundary.

        Keyword arguments:
        num_vertices -- the number of vertices the generated polygon
                        should have

    Return a list of integer tuples (pairs), its length should
    be num_vertices.
    """
    if num_vertices < 3:
        return []

    polygon = [(0,0)]
    E_start_x = 0
    current_y = -1
    current_x = 0
    add_to_x = 1

    """ adds vertices to the polygon until the length of the polygon is the
        correct length. Adds a variety of vertices so that it handles the
        3 cases in monotone_triangulate(polygon)"""
    while num_vertices > len(polygon):
        if len(polygon) % 7 == 0:
            current_x -= 1
            polygon.append((current_x, current_y))
            add_to_x = 2
            current_y -= 1

        elif len(polygon) % 5 == 0:
            E_start_x -= 1
            polygon.insert(0, (E_start_x, current_y))
            current_y -= 1

        else:
            current_x = current_x + add_to_x
            polygon.append((current_x, current_y))
            current_y -= 1
            add_to_x += 1

    return polygon

def brute_force_triangulate(polygon: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """
    Compute a triangu   lation of a polygon (given as a list of
    integer pairs (tuples), which is a clockwise ordering
    of the coordinates of the vertices) in time O(n^2). Return
    the diagonals as a list of points.

        Keyword arguments:
        polygon -- a list of integer tuples (pairs) representing
                   a clockwise ordering of points around a polygon

    Return a list of integer tuples (pairs). Supposing the list is
    called diagonals, diagonal i is stored as endpoints in
    diagonal[2i] and diagonal[2i+1].
    """
    if polygon is None or len(polygon) < 3:
        print("Error")
        return None

    if len(polygon) == 3:
        return []

    #print("current polygon")
    #print(polygon)
    v, vindex = leftmost(polygon)
    print("leftmost vertex")
    print(v)

    #reordering list so it is easier to work with
    polygon = polygon[vindex:] + polygon[:vindex]
    vindex = 0
    uindex = len(polygon) - 1
    windex = vindex + 1

    u = polygon[uindex]
    w = polygon[windex]
    print("u")
    print(u)
    print("w")
    print(w)
    z, zindex = findz(polygon, u, w, v)
    print("ZZZZZZZZZZZZZ")
    print(z)

    if z == None:
        print("v u w")
        print(v)
        print(u)
        print(w)
        print("windex uindex")
        print(windex)
        print(uindex)
        print("new shape to triangulate")
        print(polygon[windex:uindex + 1])
        return brute_force_triangulate(polygon[windex:uindex + 1]) + [u, w]


    print("two shapes to triangulate :")
    print(polygon[vindex: zindex + 1])
    print(polygon[zindex:] + [polygon[vindex]])
    return brute_force_triangulate(polygon[vindex: zindex + 1]) + [v,z] + \
    brute_force_triangulate(polygon[zindex:] + [polygon[vindex]])

def monotone_triangulate(polygon: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """
    Compute a triangulation of a polygon (given as a list of
    integer pairs (tuples), which is a clockwise ordering
    of the coordinates of the vertices) in time O(n), assuming
    that the polygon is y-monotone. Return the diagonals as a
    list of points.

        Keyword arguments:
        polygon -- a list of integer tuples (pairs) representing
                   a clockwise ordering of points around a polygon

    Return a list of integer tuples (pairs). Supposing the list is
    called diagonals, diagonal i is stored as endpoints in
    diagonal[2i] and diagonal[2i+1].
    """
    stack = [] #stack will store the indexes of the vertices in polygon list
    diagonals = []
    ordered_vertices_indexes, top_index, bottom_index = order_vertices(polygon)

    #adds the first two vertices to the stack
    stack.append(ordered_vertices_indexes[0])
    stack.append(ordered_vertices_indexes[1])

    for index in range(2, len(polygon) - 1):
        vj_index = ordered_vertices_indexes[index]
        print("current vj")
        print(polygon[vj_index])
        print("start stack")
        for i in stack:
            print(polygon[i])
        print("end stack")

        if chain_side(stack[-1], top_index, bottom_index) == \
           chain_side(vj_index, top_index, bottom_index):

            while len(stack) > 1 and not reflex(polygon, vj_index, stack, top_index, bottom_index):
                print("convex")
                print("adding diagonal")
                diagonals.append(polygon[vj_index])
                diagonals.append(polygon[stack[-2]])
                stack.pop()
            stack.append(vj_index)

        else: # v_j is lower endpoint of e
            rk_index = stack[-1]
            print("v_j is a lower endpoint of e")

            while len(stack) != 1:
                print("add diagonal")
                diagonals.append(polygon[vj_index])
                diagonals.append(polygon[stack[-1]])
                stack.pop()
            stack.pop() # toss r_1
            stack.append(rk_index)
            stack.append(vj_index)


    #processing the last vertice
    stack.pop()
    vn = polygon[ordered_vertices_indexes[-1]]
    while len(stack) != 1:
        print("end diagonals")
        diagonals.append(vn)
        diagonals.append(polygon[stack[-1]])
        stack.pop()

    return diagonals


def is_triangulation(polygon: List[Tuple[int,int]], diagonals: List[Tuple[int,int]]) -> bool:
    """
    Return True if and only if diagonals is a triangulation
    of polygon. Must be done efficiently.

        Keyword arguments:
        polygon   -- a list of integer tuples (pairs) representing
                     a clockwise ordering of points around a polygon

        diagonals -- a list of integer tuples (pairs) representing
                     line segments between vertices of the polygon
                     Diagonal i is stored as endpoints in
                     diagonal[2i] and diagonal[2i+1].

    Return True if diagonals is a trianglulation of polygon, False
    otherwise.
    """

    #must be n - 3 diagonals for there to be n - 2 triangles
    if len(diagonals) / 2 !=  len(polygon) - 3:
        return False

    #checking if any of the diagonals intersect each other
    for index in range(0, len(diagonals), 2):
        for index1 in range(0, len(diagonals), 2):
            if index != index1:
                #print(diagonals[index], diagonals[index + 1])
                #print(diagonals[index1], diagonals[index1 + 1])
                if intersection_found((diagonals[index], diagonals[index + 1]), \
                (diagonals[index1], diagonals[index1 + 1])):
                    print((diagonals[index], diagonals[index + 1]))
                    print((diagonals[index1], diagonals[index1 + 1]))
                    print("HERE")
                    return False

    return True

def test() -> None:
    """ Test your triangulation routines. """
    #print(brute_force_triangulate([(4,0), (2,2), (5,4), (6,2), (3,2.5)]))
    #print(brute_force_triangulate([(-6,0), (-8,0), (-10,4), (-9,6), (-7,4), (-4,2)]))
    #print(brute_force_triangulate([(-7, 4), (-4, 2), (-6, 0), (-8, 0), (-9, 6)]))
    #print(brute_force_triangulate([(2,2), (8,6), (4,3), (6,3), (4,2), (9,-1), (10,-2)]))
    #print(brute_force_triangulate([(2,0), (3, 3), (4,4), (1,-2)]))
    #print(is_triangulation([(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)], [(1,4), (1,5), (-4,3),(2,2), (5,5), (6,6)]))
    #issue polygon
    #polygon = [(-1,-3), (-2, -2), (0,0), (2,-3), (5, -5), (-2,-15), (-6,-6)]
    #polygon = [(-2,-5), (0,0), (1,-2), (2.5,-3), (4, -3.5), (5,-3.7), (-2,-5)]
    #print(monotone_triangulate(polygon))    #print(brute_force_triangulate(generate_monotone_polygon(20)))

    #polygon = generate_monotone_polygon(10)
    #polygon = [(0,0), (1, -2), (3,-3), (2.5, -4), (4, -5), (5, -7), (0, -10), (-3, -8), (-2, -6)]
    #print(polygon)
    polygon = [(-1,-5),(0,0),(1,-1),(3,-2),(6,-3),(10,-4),(15,-6),(14,-7),(16,-8),(19,-9)]
    #diagonal1 = brute_force_triangulate(polygon)
    diagonal1 = monotone_triangulate(polygon)
    print(diagonal1)
    if not is_triangulation(polygon, diagonal1):
        print("error is monotone_triangulate")

    #if not is_triangulation(polygon, diagonal2):
    #    print("error in monotone_triangulate")


if __name__ == '__main__':
    test()
