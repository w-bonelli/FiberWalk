# FiberWalk

The Fiberwalk demo depends on:
- Python 2.7
- [`networkx`](http://networkx.lanl.gov/) 2.2
- [`numpy`](http://sourceforge.net/projects/numpy/) 1.16
- [`scipy`](http://www.scipy.org/SciPy) 1.2.3
- [`matplotlib`](https://matplotlib.org/) 2.2.5
- the modules `lattice` and `fiberWalk`

Older library versions may work, but are not guaranteed to.

The graphics shown in the paper where generated with [mayavi2](http://docs.enthought.com/mayavi/mayavi/index.html).

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the FiberWalk paper if you use the code for your project.
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0085585

# Author

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


## Running the demo

1.) make sure you have python 2.7, numpy 1.3, scipy 0.9 and networkX 1.7 properly installed
The Demo may run with older versions, but no guarantee 

2.) copy the files Lattice.py, FiberWalk.py and main.py to a folder your choice.

3.) on the command line, go to the folder of your choice

4.) type: python main.py

5.) The program generates saves plots in the creates subfolder FiberWalks and python pickles of the walk data in the folder you have chosen.

6.) Happy Walking!

Changing walk parameters:

2.) search for the following section

### Configure Parameters    

Open main.py in a text editor, then search for the following variables:

`dimension`: Number of dimensions, defaults to 2.
`numberOfSteps`: Length of the walk, defaults to 100 steps.
`numberOfObjects`: Number of walks to simulate, defaults to 10.
`numberOfContractions`: Number of contractions per step, defaults to 1.

3.) adjust the variables as you wish

4.) save the file and run the demo

5.) Happy Walking!
