# Adding Dynamic Vertices to a Dynamic de Bruijn graph.

Belazzougui et al. (2016) presents  a space- and time-efficient fully dynamic implementation de Bruijn graphs. The goal of this project is to implement the approach in this paper for dynamic vertices. 

## Getting Started

You can find the project code hosted as a zip file on the project website: https://h-bee.github.io/ or you can choose to clone it from the GitHub repository : https://github.com/thegreedychoice/FullyDynamicDeBruijn

### Prerequisites

You will need:

* Python 2.7 (2.7.12 or newer)
* graphviz (an open source visualization software)
* BitArray2D (for storing the IN and OUT matrix efficiently)

### Installing

* Download and install **Python version (2.7.12 or newer)** from [Python's official download page] (https://www.python.org/downloads/) if it's not present.
* Install **pip** using this [link](https://pip.pypa.io/en/stable/installing/) if it's not already installed.
* Install **graphviz** using the following command:
```
    
                pip install graphviz

```
* Install **BitArray2D** using the following command:
```

                pip install BitArray2D

```
## Running the program

For testing purposes we have created 3 instances of varying file sizes(around 5MB, 10MB, 20MB) FASTQ file (s_6_1.fastq.gz) obtained from [here] (http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_mc/). 

Run **'test.py'** and follow prompts to input size of K and k-mers for searching, insertion and deletion.
```

            python test.py
```


### Changing the input FASTQ file

### Output Specifications




## Built With

* [Python](https://www.python.org/doc/) - The programming language used
* [graphviz](https://www.graphviz.org/documentation/) - Visualization package for graph
* [BitArray2D](https://pypi.python.org/pypi/BitArray2D/2.1) - Used for memory-efficient packed representation of 2D bit arrays in Python
* [Minimal Perfect Hashing] (http://stevehanov.ca/blog/index.php?id=119) - Used for performing minimal perfect hashing on top of the Rabin-Karp hashed values

## Versioning

We use [GitHub](https://github.com/) for versioning. For the current version, see the [Final project code for submission](https://github.com/thegreedychoice/FullyDynamicDeBruijn/tree/v1.0). 

## Authors

* **Anthony M Colas**
* **Ayush Khandelwal** 
* **Harish Balaji**
* **Shubham Shukla**

## Acknowledgments

* **Dr. Christina Boucher**, for guiding us with the project all throughout the course.
* **Dr. Travis Gagie**, for providing insights on his paper 'Fully dynamic de Bruijn Graph' and clearing doubts regarding forest data structure used  
* **Steve Hanov**, for the code for an easy, minimal perfect hashing method hosted on his blog. 

