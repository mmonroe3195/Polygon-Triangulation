"""
    File:        triangulate_timer.py
    Author:      Darren Strash
    Course:      CS 307 - Computational Geometry
    Assignment:  Problem Set 4 - Polygon Triangulation
    Description: Execute different triangulation algorithms (implemented by
    you!) and print out a table of times to show growth rates.
"""

from triangulate import generate_monotone_polygon
from triangulate import brute_force_triangulate
from triangulate import monotone_triangulate

import datetime

NUM_ITERATIONS = 3

def average_time(algorithm, polygon):
    """call algorithm on polygon repeatedly and return average time"""
    time_diff = None

    for iteration in range(0, NUM_ITERATIONS):
        polygon_copy = polygon[:]
        start = datetime.datetime.now()
        algorithm(polygon_copy)
        end = datetime.datetime.now()
        if time_diff == None:
            time_diff = end - start
        else:
            time_diff = time_diff + end - start

    time_ms = time_diff.microseconds // 1000
    time_ms = time_ms + time_diff.seconds * 1000
    time_ms = time_ms + time_diff.days * 24 * 60 * 60 * 1000

    return time_ms

def get_table_entry(num_vertices, item):
    """get the appropriate table entry, which is either a number of vertices
       or a running time"""
    polygon = generate_monotone_polygon(num_vertices)
    if item == "n":
        return num_vertices
    elif item == "brute force":
        return average_time(brute_force_triangulate, polygon) 
    elif item == "monotone":
        return average_time(monotone_triangulate, polygon) 
        
    return -1

def build_header_and_legend(option):
    """construct the header entries, which are also used to fill table entries"""
    # always print n (number of vertices) and running time of monotone algorithm
    header = ["n", "monotone"]

    print("Legend:")
    print("  n           : the number of vertices in the polygon")
    print("  monotone    : the running time of the monotone triangulation algorithm (in ms)")

    if option == "all":
        header.append("brute force")
        print("  brute force : the running time of Brute Force triangulation (in ms)")

    print("")

    return header

def run_experiment(option):
    """run the timing experiement according to the user-supplied option"""
    header = build_header_and_legend(option)

    for item in header:
        print("{:>15} ".format(item), end="")
    print("")

    for i in range(2,35):
        size = 2**i
        for item in header:
            print("{:>15} ".format(get_table_entry(size, item)), end="")
        print("")

def main():
    """Get user input and run appropriate timing experiment."""
    print("Welcome to Triangulate Timer! Press Ctrl+C at any time to end...")

    option = input("Which test would you like to run (all,monotone)? ")
    while option not in ["all", "monotone"]:
        print("Unrecognized option '", option, "'")
        option = input("Which test would you like to run (all,monotone)? ")

    if option == "all":
        print("Running all algorithms with n vertices.")
        print("This test includes the brute force algorithm. Run 'monotone' to remove it.")
    else:
        print("Running all algorithms except the bruteforce algorithm, with n vertices.")

    run_experiment(option)

if __name__ == "__main__":
    main()
