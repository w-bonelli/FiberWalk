# FiberWalk

The Fiberwalk demo depends on:
- Python 2.7
- [`networkx`](http://networkx.lanl.gov/) 2.2
- [`numpy`](http://sourceforge.net/projects/numpy/) 1.16
- [`scipy`](http://www.scipy.org/SciPy) 1.2.3
- [`matplotlib`](https://matplotlib.org/) 2.2.5
- the FiberWalk modules `lattice` and `walk`

Older library versions may work, but are not guaranteed to.

The graphics shown in the paper where generated with [mayavi2](http://docs.enthought.com/mayavi/mayavi/index.html).

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

1) Make sure you have Python 2.7, `numpy` 1.16, `scipy` 1.2.3, `matplotlib` 2.2.5, and `networkx` 2.2 properly installed. The Demo may run with older versions, but no guarantee.

2) From the project root, run `python main.py`. Happy walking!

## Author

Author: Alexander Bucksch
Department of Plant Biology
Warnell School of Forestry and Natural Resources
Institute of Bioinformatics
University of Georgia, Athens

Mail: bucksch@uga.edu
Web: http://www.computational-plant-science.org

## License

Copyright (c) 2012, 2016 Alexander Bucksch
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the Fiber Walk Demo Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

'''
