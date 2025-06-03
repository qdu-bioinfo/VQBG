# VQBG

VQBG is a powerful tool for virus genome assembly. The package includes precompiled binaries for direct execution.

## Requirements

- **Operating System:** Linux
- **C++ Standard:** C++11 or higher
- **Boost**
```
 wget http://downloads.sourceforge.net/project/boost/boost/1.80.0/boost_1_80_0.tar.gz
 tar xfz boost_1_80_0.tar.gz
 rm boost_1_80_0.tar.gz
 cd boost_1_80_0
 ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options,regex,filesystem,system
 export
 ./b2 install
 cd /home
 rm -rf boost_1_80_0
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
```
## Installation
### 1. Install via Bioconda
```
conda install vqbg
```
### 2.Use Precompiled Binaries:
Download VQBG from GitHub or clone the repository:
```
git clone https://github.com/qdu-bioinfo/VQBG.git
```
 Navigate to the FC-Virus directory
 ```
cd /path/to/VQBG
```
 Install by running make
 ```
make
```
 Set up environment variables
 ```
vim ~/.bashrc
```
 Add the following line to ~/.bashrc
 ```
export PATH="$PATH:/path/to/VQBG"
```
 Apply the changes
 ```
source ~/.bashrc
```
### 3.Compile from source
To compile from source, navigate to the code directory:
```
cd /path/to/VQBG
```
Compile using g++:
```
g++ -o VQBG common.cpp kmer_hash.cpp main.cpp sequence_graph.cpp utility.cpp
```
## Running VQBG
VQBG supports assembly results from [FC-Virus](https://github.com/qdu-bioinfo/fc-virus).
### Quick Usage
```
===============================================================================
Usage: Assemble [--reads/--kmers] <filename>  [opts]
===============================================================================
**Required :
--reads/-i <string>           : the name of the file containing reads

** Optional :
--kmer_length/-k <int>        : length of kmer, default: 25.
--fasta/-a                    : input reads file is in fasta format.
--fastq/-q                    : input reads file is in fastq format.
--output_filename/-o <string> : Name of the output file, default: paths.fasta.
--trunk_filename/-t <string> : Name of the trunk file.
--help/-h                     : display the help information.

================================================================================

```
### Then run VQBG:
```
./VQBG -k 25 -q -t FC-Virus.fa -o VQBG.fasta -i forward.fastq -i reverse.fastq
```

### Example command
```
./VQBG -k 25 -q -t ./path/to/FC-Virus.fa -o VQBG.fasta -i ./path/to/forward.fastq -i ./path/to/reverse.fastq
```
## Experiment

1.Simulated Dataset, can be found at [savage-benchmark](https://bitbucket.org/jbaaijens/savage-benchmarks/src/master/) 

- 6 Poliovirus (20,000x)

- 10 HCV (20,000x)

- 15 ZIKV (20,000x)

2.Simulation data with different sequencing depth and different strain abundances, can be found at [VQBG-data](https://bitbucket.org/vqbg-benchmark/data/src/main/)

- 6 Poliovirus (6,000x)

- 6 Poliovirus (200x-14,000x)

3.Real Dataset

- 5 HIV labmix (20,000x) [SRR961514](https://www.ncbi.nlm.nih.gov/sra/?term=SRR961514), reference genome sequences are available at [5 HIV References](https://github.com/cbg-ethz/5-virus-mix/blob/master/data/REF.fasta)

- 2 SARS-COV-2 (4,000x) [SRR18009684](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009684), [SRR18009686](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18009686), reference genome sequences are available at [2 SARS-COV-2 Dataset](https://github.com/RunpengLuo/sarscov2-4000x)
