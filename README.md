# Problem Set 4: Polygon Triangulation

## Setup

1. From the terminal, clone your repository with the command

    `git clone your_repository_link`

    where `repository_link` is likely of the form ```https://github.com/HamiltonCollege/ps04-triangulate-username```.

2. `cd` into the directory created (of the form ps04-triangulate-*username*), and replace or update *triangulate.py* with your functions from the assignment.

3. From the repository directory, run the timing code with 

    `python3 triangulate_timer.py`

    then type in one of the two experiment options

    `all monotone`

    which will run a subset of the algorithms so that you can evaluate the running time.

## Evaluating growth rates

**Experiment**: `all`

In this experiment a single table is printed, one line at a time. On each line of the experimental output, the number of vertices *n* on the polygon doubles, and the running times for the monotone and brute force triangulation algorithms are printed. Evaluate the growth rate of your algorithms by looking at how the running time grows with the number of vertices:

- For *O(n)* running time, you should see the running time roughly double as *n* doubles.
- For *O(n^2)* running time, you should see the running time roughly quadruple as *n* doubles.

**Experiment**: `monotone`

In this experiment, the brute force algorithm is not included. Use this experiment to evaluate growth rates of the y-monotone triangulation algorithm for sizes of polygons that the brute force algorithm cannot handle in a reasonable amount of time.

## Submitting

You only need to submit the file *triangulate.py*. The easiest way is to upload it to GitHub by opening your repository in a browser (you can find it by navigating to github.com and clicking the link on the left-hand side); then click the *Upload files* button, and drag *triangulate.py* to the browser window to update it.

Otherwise, from your repository directory on your machine, use `git` commands to submit your code with: 

`git add triangulate.py`

`git commit -m "commit message goes here"`

`git push origin master`
