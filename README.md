# FiberWalk

The Fiberwalk demo depends on:
- Python 3.7
- [`networkx`](http://networkx.lanl.gov/) 2.4
- [`numpy`](http://sourceforge.net/projects/numpy/) 1.18.4
- [`mayavi`](https://docs.enthought.com/mayavi/mayavi/mlab.html) 4.7.1
- `PyQt5` 5.14.2
- the FiberWalk modules `lattice` and `walk`

Older library versions may work, but are not guaranteed to.

The graphics shown in the paper where generated with [`mayavi`](http://docs.enthought.com/mayavi/mayavi/index.html).

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the FiberWalk paper if you use the code for your project.
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0085585

## Running the demo

### Parameters    

Open `demo.py` in a text editor, then search for the following variables:

- `dimensions`: Number of dimensions, defaults to 2.
- `steps`: Length of the walk, defaults to 100 steps.
- `walks`: Number of walks to simulate, defaults to 10.
- `contractions`: Number of contractions per step, defaults to 1.

Then save the file.

### Run

1) Make sure you have Python 2.7, `numpy` 1.16, `scipy` 1.2.3, `matplotlib` 2.2.5, and `networkx` 2.2 properly installed. The demo may run with older versions, but no guarantee.

2) From the project root, run `python demo.py`. Happy walking!

## Author

Alexander Bucksch

Department of Plant Biology

Warnell School of Forestry and Natural Resources

Institute of Bioinformatics

University of Georgia, Athens

## Contact

- Email: bucksch@uga.edu
- Web: http://www.computational-plant-science.org
