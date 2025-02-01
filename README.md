# Chemical Equitable Partitions: Algorithm and Implementation
 
*Implementation in Fortran, of a well-known algorithm that refines an existing 
graph partition until finding the coarsest equitable partition. The algorithm is
adapted to the problem of finding the 'Chemical' Equitable Partition of a given 
molecular graph, i.e. the coarsest refinement of an initial partition based on 
the chemical elements and valences of the atoms or atom groups of the molecule.*  

#### Author:

Vasilios Raptis  
  
#### Contact: 

    v.raptis@external.euc.ac.cy  
  
**Version**        :   1.0.0  

**Release date**   :   January 22, 2025  


## 1. Introduction
This project contains part of the software that was developed and employed to 
calculate the results reported in the article "Graph Partitions in Chemistry"
by I. Michos and V. Raptis (Entropy, 25, p. 1504, 2023, doi: 10.3390/e25111504)
The project consists of a program in Fortran, which reads an input file written 
in DOT, and describing a molecular graph; an auxiliary Fortran module borrowed 
from another project of mine; a simple shell script to automate menial tasks; 
examples of input files; and this README file.  

In the below paragraphs, a Linux environment (bash shell) is assumed. It should 
be very easy to adapt the instructions to other environments.  

## 2. How to build and run the program 

### Compilation 
The source code is in Fortran so you need the appropriate compiler. GNU Fortran 
is generally available for installation on most common Linux distributions. For 
example on Debian or Ubuntu systems, you can install it (if not already there) 
with the following commands:  
    
    sudo apt-get update
    sudo apt-get install gfortran
    
Once you have GNU Fortran in place, compilation is straightforward.  

    gfortran -c strings.f90
    gfortran -c cep.f90
    gfortran strings.o cep.o -o cep.exe

The above commands will generate an executable called *cep.exe* in your folder.  

### Running the program 
Assuming an input file called *example.dot* in place, you just have to run  
 
    ./cep.exe example.dot

and you should see something like the following on your screen:  

     ncell= 2
     ncell= 3
     ncell= 4
     ncell= 4
    
    *** Chemical Equitable Partition found! ***
        Number of cells                 :    4
        % compression ratio             :   50.00
        Information content             :    2.00
        % normalised information content:   66.67
        
(The number of lines of the form *ncell=...* and the exact values printed on 
screen, will vary with the graph described in the input file.)  

With the end of the calculation, two more files are generated, called *qg.dot* 
and *qg.dat*. These are discussed in more detail in Section 4.  

The sought partition and its quotient graph come about as the refinement of an 
initial partition based on the chemical elements and valences of atoms in the 
original molecular graph. In the improbable case of failure of the refinement 
algorithm to converge after 1000 trials, execution will stop and an informative 
message will appear on screen.  

## 3. Input explained 
Input is read from files written in DOT, a very simple and intuitive language 
for the description and visualisation of graphs. DOT is commonly employed with 
the GraphViz open source toolkit. More information about DOT and GraphViz, can 
be found on the latter's website: https://graphviz.org  

## 4. Output explained 

### File 'qg.dot' 

This file is written in DOT and describes the quotient graph associated with the 
graph's chemical equitable partition. Once we are familiar with DOT, reading and
understanding qg.dot is straightforward.  

### File 'qg.dat'  

This file summarises the information pertinent to the graph's chemical equitable 
partition, namely the converged connectivity profile of each graph node; lists 
of the nodes contained in each partition cell; the quotient matrix; and lastly, 
a summary of the metrics pertaining to the partition operation (number of cells, 
compression ratio, graph information content and its normalised counterpart).  

## 5. The companion script 

This script automates the procedure of calling the appropriate GraphViz utility 
to read a .dot file and visualise the corresponding graph by generating images 
in png and svg format.  

Assuming an input file called *example.dot* in place, you just have to run  
 
    ./graph.sh example.dot

and the files *example.png* and *example.svg* will be generated.  

Internally, the script reads the appropriate keyword inside the .dot file and 
determines whether to call the *neato* or *dot* utility (suitable for undirected 
or directed graphs, respectively).  

## 6. A glance at the program's internals  
The program implements a simple algorithm (described in reference 1, section 7)
which has long been known as an intermediate step of existing algorithms built 
for the problem of graph isomorphism - see more references in section 7.  

The first step is to read the dot file, save nodes and edges listed therein, and 
form the adjacency matrix.  

Then, the connectivity profiles of the nodes are initialised. These are arrays 
that are readjusted dynamically to accomodate the number of neighbours each node 
has in a given partition cell. Therefore, the number of array elements changes 
to account for more cells as the partition is being refined.  

At the beginning, connectivity profiles contain only node degrees and chemical 
identities because the initial partition to refine is based on that information, 
but more array elements are to be added gradually.  

Then, follows a loop that stops when the partition cells defined by the varying 
connectivity profiles with each step, converge to a stable configuration. With 
each new step, the old partition is saved; a new refinement is defined such that
each cell contains nodes with identical connectivity profiles; and the refined 
partition is compared to the old one. When the coarsest equitable partition that 
refines the initial partition, is found, the loop is terminated.  

When the algorithm converges, the quotient graph is determined by populating its
adjacency matrix and writing the corresponding dot file.  

Finally, metrics describing the operation are calculated and printed onscreen as 
well as in a dedicated file, namely: 
* number of converged partition cells 
* compression ratio 
* information content 
* normalised information content 

## 7. Literature

1.  I. Michos, V. Raptis. Graph Partitions in Chemistry. *Entropy*, **2023**, *25* 
    (11), 1504.  
2.  Neave O’ Clery, Ye Yuan, Guy-Bart Stan, Mauricio Barahona. Observability and 
    coarse graining of consensus dynamics through the external equitable partition. 
    *Phys. Rev. E* **2013**, *88*, 042805.  
3.  M. Cao, S. Zhang, and M. Camlibel. A Class of Uncontrollable Diffusively 
    Coupled Multiagent Systems with Multichain Topologies. *IEEE Transactions on 
    Automatic Control* **2013**, *58*, 465.  
4.  Sergey Krivovichev. Topological complexity of crystal structures: quantitative
    approach, *Acta Cryst.* **2012**, *A68*, 393.
5.  M. Dehmer, A. Mowshowitz. A history of graph entropy measures. *Inf. Sci.*,
    **2011**, *181*, 57-78.
6.  M. Dehmer, A. Mowshowitz. Inequalities for entropy-based measures of network
    information content. *Appl. Math. Comput.*, **2010**, *215*, 4263–4271. 
7.  G. Hahn and G. Sabidussi. Graph Symmetry: Algebraic Methods and Applications, 
    NATO ASI Series. Series C, Mathematical and Physical Sciences, Vol. 497 
    (Kluwer Academic Publishers, Dordrecht, 1997).  
8.  Brendan D. McKay. Practical Graph Isomorphism. *Congressus Numerantium*, **1981**, 
    *30*, 45–87.  
9.  D. Cvetkovic, M. Doob, H. Sachs.  Spectra of Graphs: Theory & Application 
    (Volume 87 of Pure and Applied Mathematics). Academic Press, 1980.  
10. Brendan D. McKay. Computing automorphisms and canonical labellings of graphs. 
    *Combinator. Math.* **1978**, *686*, 223–232.  
11. D. G. Corneil and C. C. Cotlieb. An efficient algorithm for graph isomorphism, 
    *J. Assoc. Comput. Mach.* **1970**, *17*, 51–64.
12. A. Mowshowitz. Entropy and the Complexity of Graphs: IV. Entropy Measures and
    Graphical Structure. *Bull. Math. Biophys.*, **1968**, *30*, 533-546.
13. A. Mowshowitz. Entropy and the Complexity of Graphs: III. Graphs with Prescribed
    Information Content. *Bull. Math. Biophys.*, **1968**, *30*, 387-414.     
14. A. Mowshowitz. Entropy and the Complexity of Graphs: II. The Information Content
    of Digraphs and Infinite Graphs. *Bull. Math. Biophys.*, **1968**, *30*, 225-240. 
15. A. Mowshowitz. Entropy and the Complexity of Graphs: I. An Index of the Relative
    Complexity of a Graph. *Bull. Math. Biophys.*, **1968**, *30*, 175-204. 
16. E. Trucco. A Note on the Information Content of Graphs. *Bull. Math. Biophys.*,
    **1956**, *18*, 129-135. 
17. N. Rashevsky. Life, Information Theory, and Topology. *Bull. Math. Biophys.*,
    **1955**, *17*, 229-235. 

## 8. Legal stuff, etc.
This program implements an established algorithm to find the coarsest equitable 
partition that refines an existing partition of a given graph. It is freely 
distributed under the MIT license.  

### How to cite 
Please cite as follows:  

Vasilios Raptis, 
Chemical Equitable Partitions: Algorithm and Implementation, 
https://github.com/brapt2005/cep (accessed on: ... insert date...)  

You may also cite the following article:  

Michos I, Raptis V  
'Graph Partitions in Chemistry', Entropy, **2023**, *25* (11), 1504  
DOI: 10.3390/e25111504  
URL: https://www.mdpi.com/1099-4300/25/11/1504  


### License

Copyright (c) 2025 Vasilios Raptis <v.raptis@external.euc.ac.cy>  

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.  

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.  

