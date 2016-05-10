# KCMBT: A very fast _k_-mer counter

###### Copyright 2015 Abdullah-Al Mamun <br />
abdullah.am.cs (at) engr.uconn.edu

 	KCMBT is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    KCMBT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KCMBT. If not, see <http://www.gnu.org/licenses/>.

## What is KCMBT
KCMBT (_k_-mer Counter based on Multiple Burst Trees) is a very fast multi-threaded _k_-mer counting algorithm. It uses cache efficient burst tries to store _k_-mers. Experimental results show that it outperforms all well-known algorithms.

## Compilation
To compile the source code, type

```
make
``` 

It will create three executable files in bin directory:

* _kcmbt_ generates binary files containing _k_-mers with their counts
* _kcmbt\_dump_ produces human readable file having _k_-mers with their counts
* _kcmbt_query_ copies range specific _k_-mer list to user-defined output file

## Usage
Please run ./bin/kcmbt at first, and then run ./bin/kcmbt\_dump. ./bin/kcmbt\_dump uses kcmbt generated files as input, and outputs list of human readable _k_-mers with their counts.

To run kcmbt, use

```
./bin/kcmbt -k <k-mer length> -i <@file_listing_fastq_files or fastq_file> -t <number_of_threads>
Parameters:
	-k:	k-mer length (10 <= k <= 32, default 28) 
	-i:	input file in fastq format (start with @ if the file contains a list of fastq files)
	-t:	number of threads (please use 2^x threads, x = 0, , 2, 3, ..)
```
Example: ```./bin/kcmbt -k 28 -i srr.fastq -t 4```

To run kcmbt_dump, use

```
./bin/kcmbt_dump number_of_threads_used_in_kcmbt
```
Example: ```./bin/kcmbt_dump 4```

_kcmbt\_dump_ creates kmer\_list.txt which contains human readable _k_-mer list with their counts. This file may be huge. So we can query a specific range of _k_-mer lists using _kcmbt\_query_.

To run kcmbt\_query, use

```
./bin/kcmbt_query out_file_name begin_count end_count
```
Example: ```./bin/kcmbt_query out 2 50```

## Contact
For questions, suggestions, bugs, and other issues, please contact:

```
Abdullah-Al Mamun
abdullah.am.cs (at) engr.uconn.edu
```
