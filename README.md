This project is a fork of <a href="https://github.com/karmaresearch/vlog">VLog</a>.

[![Build Status](https://travis-ci.org/karmaresearch/vlog.svg?branch=master)](https://travis-ci.org/karmaresearch/vlog)

## Installation 

We used CMake to ease the installation process. To build GLog, the following
commands should suffice:

```
mkdir build
cd build
cmake ..
make
```

External libraries should be automatically downloaded and installed in the same
directory. The only library that should be already installed is zlib, which is
necessary to read gzip files. This library is usually already present by
default.

To enable the web-interface, you need to use the -DWEBINTERFACE=1 option to
cmake.

If you want to build the DEBUG version of the program, including the web
interface: proceed as follows:

```
mkdir build_debug
cd build_debug
cmake -DWEBINTERFACE=1 -DCMAKE_BUILD_TYPE=Debug ..
make
```

## VLDB 2021

To facilitate the reproduction of the experiments presented in the paper

> Tsamoura, Efthymia, David Carral, Enrico Malizia, and Jacopo Urbani. "Materializing knowledge bases via trigger graphs." Proceedings of the VLDB Endowment 14, no. 6 (2021): 943-956

we have copied all the datasets, scripts, and other useful resources in the folder:

https://drive.google.com/drive/folders/14ZIn8fWkZ7oCbdGZXBTUI7kiI-1TD_kM?usp=sharing

Below, re report more detailed instructions. Please be aware that the
instructions refer to the version of the code used at the time of the paper. If
you download the last version of GLog, it might be that the syntax (or the
performance) has changed.

## Instructions to replicate the experiments

### Linear scenarios ###

The folder "Linear " includes the scripts and the KBs to reproduce the
experiments for the linear scenarios. The KBs for each scenario are included
in the folders named **lubm125**, **uobm**, **dbpedia**, **claros** and
**Reactome**.

Each folder includes:
* a **-linear.dlog** file that includes the linear rules of the KB;
* a folder named **db** that includes the extensional data in .csv format;
* a filed named **edb.conf** that specifies the path to the extensional data.

#### Running GLog for the linear scenarios 

To materialize a KB using VLog run the command:

```
./vlog mat --edb scenario/edb.conf --rules scenario/program-linear.dlog --storemat_path derived_data_folder 2>&1 | tee logs.txt
```

* **scenario** is the top-level folder of each scenario;
* **edb.conf** and **program.dlog** are as above;
* **derived_data_folder** is the folder to store the derived data; and
* **logs_folder/log** is the file to store the log file.

#### Running TG-guided materialization for the linear scenarios 

To materialize a KB using TGs run the commands below:

```
./vlog trigger --edb scenario/edb.conf --rules scenario/program-linear.dlog --trigger_paths tgfile --trigger_algo linear 2>&1 | tee logs.txt
```

The above command computes and minimizes a TG using the _tglinear_ and the _minLinear_ algorithms from Section 5 (option **--trigger_algo linear**)
and stores it in the file named _tgfile_ (option **--trigger_paths tgfile**).

```
./vlog  mat_tg --edb scenario/edb.conf --rules scenario/program-linear.dlog --trigger_paths tgfile --storemat_path derived_data_folder 2>&1 | tee logs.txt
```

The above command consumes the previously computed TG (option
**---trigger_paths tgfile**) and performs TG-guided materialization.

#### Examples for the linear scenarios

The script **launch_experiments.sh** includes the commands to reproduce the
experiments for the linear scenarios.

**Notice:**  Launch the commands from the directories where the **edb.conf**
files are, so that the relative paths are properly parsed.

### Datalog scenarios

The folder "Datalog" includes the scripts and the KBs to reproduce the
experiments for the Datalog scenarios.  The KBs for each scenario are included
in the folders named **lubm_125**, **uobm**, **dbpedia** and **claros**.

Each folder includes:
* an **_L.dlog** file that includes the "L" rules of the KB;
* an **_LE.dlog** file that includes the "LE" rules of the KB;
* a folder named **db** that includes the extensional data;
* a filed named **edb.conf** that specifies the path to the extensional data.

#### Running VLog for the Datalog scenarios 

To materialize a KB using VLog run the command:

```
./vlog mat --edb scenario/edb.conf --rules scenario/program.dlog --storemat_path derived_data_folder --decompressmat 1 2>&1 | tee logs.txt
```

The options are as for the linear scenarios. The option **--decompressmat 1** de-compresses the derived data.

#### Running TG-guided materialization for the Datalog scenarios

To materialize a KB using TGs run the command:

```
./vlog tgchase --edb scenario/edb.conf --rules scenario/program.dlog --querycont flag1 --edbcheck flag2 --storemat_path derived_data_folder --decompressmat 1 2>&1 | tee logs.txt
```

Set **flag1** to **1** to minimize a TG using the _minDatalog_ algorithm from Section 6.1; otherwise set **flag1** to **0**.

Set **flag2** to **1** to apply the optimized rule execution strategy from Section 6.2; otherwise set **flag2** to **0**.

#### Examples for the Datalog scenarios 

The script **run.sh** inside the folder "Datalog" includes the commands to reproduce the experiments for the Datalog scenarios.

**Notice:**  Launch the commands from the directories where the **edb.conf** files are, so that the relative paths are properly parsed.

### RDFS scenarios ###

The folder "RDFS" includes the scripts and the KBs to reproduce the experiments for the RDFS scenarios.
Inside "RDFS" there are three subfolders:
**vlog_glog**;
**inferray**; and
**webpie**.

Each subfolder includes the scripts and the data to run VLog, GLog, Inferray and WebPIE.

#### Running VLog and GLog for the RDFS scenarios 

Follow the instructions for the Datalog case.
The files **lubm-edb.conf** and **rhordf_rules_from_LUBM_L** specify the path to the extensional data and the rules for the LUBM benchmark,
while the files **yago-edb.conf** and **rhordf_rules_from_YAGO** specify the path to the extensional data and the rules for the YAGO benchmark.

#### Running Infearray for the RDFS scenarios

To materialize using Inferray run the command:

```
java -Xmx32G -cp inferray-core-0.0.2-SNAPSHOT-jar-with-dependencies.jar fr.ujm.tse.lt2c.satin.inferray.benchmark.InferraySpeed lubm125_cleaned.nt 2>&1 | tee inferrary-lubm.log

java -Xmx32G -cp inferray-core-0.0.2-SNAPSHOT-jar-with-dependencies.jar fr.ujm.tse.lt2c.satin.inferray.benchmark.InferraySpeed yago-1.0.0.nt 2>&1 | tee inferrary-yago.log
```

where **lubm125_cleaned.nt** and **yago-1.0.0.nt** are the data in triples for the LUBM and YAGO benchmarks (these files are available in the folder "raw_data").


#### Running WebPIE for the RDFS scenarios 

To materialize using WebPie run the command:

```
java -Xmx32g -cp webpie-all.jar lubm125 --fragment rdfs 2>&1 webpie-lubm.log

java -Xmx32g -cp webpie-all.jar yago --fragment rdfs 2>&1 webpie-yago.log
```

where **lubm125** and **yago** are folders including the data for the LUBM and YAGO benchmarks.


### chaseBench scenarios ###

The folder "chaseBench" includes the jars, the scripts and the KBs to reproduce the experiments for the chaseBench scenarios. Each folder includes:
* a **rules.dlog** file with the rules of the KB;
* a filed named **edb.conf** that specifies the path to the extensional data;
* a folder **data** including the data in .csv format.


#### Running VLog for the chaseBench scenarios

To materialize a KB using VLog run the command:

```
./vlog mat --edb edb.conf --rules rules.dlog 2>&1 | tee logs.txt
```

#### Running TG-guided materialization for the chaseBench scenarios 

To materialize a KB using TGs run the command:

```
./vlog tgchase --edb edb.conf --rules rules.dlog --querycont 0 --edbcheck 0 2>&1 | tee logs.txt
```

#### Examples for the chaseBench scenarios 

The script **run.sh** inside the folder "chaseBench" includes the commands to reproduce the experiments for the chaseBench scenarios.


### RDFox

The folder "rdfox" includes the jars, the scripts and the KBs to reproduce the
experiments for the linear, Datalog and the two benchmarks from
**https://github.com/dbunibas/chasebench**, STB-128 and Ontology-256, using
RDFox. The benchmarks are inside the folder "scenarios". Notice that the KBs
are in the chaseBench format (described in the paper *Benchmarking the Chase*,
PODS 2017), since the RDFox jars that are available in
**https://github.com/dbunibas/chasebench/tree/master/tools/rdfox** require
inputs in this format.

Inside the folder "rdfox/scenarios" there are the subfolders "chaseBench",
"Linear" and "Datalog". These subfolders include the KBs of the corresponding
scenarios. Each scenario e.g., Claros-L, has a folder including the extensional
data in a single zip file named "data.zip". The contents of the zip files must
be extracted into the data folder before running the experiments.

### Converting the scenarios to the chaseBench format 

chaseBench requires the data in the following format:
* a file for the schema of the source relations;
* a file for the schema of the target relations;
* a file with rules from the source relations to target relations;
* a file with rules involving only taget relations;
* a folder with the data of the source schema in csv format.

To convert the triples in this format, we followed the steps below:
* we used the file **nt2csv.py** (also available in the repo) to convert the triples to the corresponding set of rules and database data;
* we created an empty file for the source schema;
* we created an empty file for the source to target rules;
* we created a target schema file and included the definitions of all relations appearing in all rules extracted in the first step;
* we created a target rules file includes all rules extracted in the first step.

Notice that the data inside the folder "rdfox/scenarios" is already in the right format and no conversion is required.

### Running RDFox 

To materialize a KB in the chaseBench format using RDFox run the command (in a Linux or Mac platform):

```
java -jar chaseRDFox.jar \
-threads <number of threads> \
-chase skolem   \
-s-sch   <path to the source schema .s-schema.txt>  \
-t-sch   <path to the target schema .t-schema.txt>  \
-st-tgds <path to the source to target dependencies .st-tgds.txt>  \
-t-tgds  <path to the target dependencies .t-tgds.txt>  \
-src     <path to the source to the data folder>
```

#### Examples for running RDFox

The three scripts **RDFox-.sh** inside the folder "rdfox" include the commands to reproduce the experiments using RDFox.

#### Running TG-guided materialization without the optimization from [1] 

To run TG-guided materialization without the rewriting technique from [1], use the argument

```
--rewritecliques 0
```

Notice that by default GLog employs this optimization.


### Scalability experiments

To reproduce the experiments on the scalability, first create the LUBM KBs using the official generator from **http://swat.cse.lehigh.edu/projects/lubm/**. With the generator, we can create a KB with a variable number of universities. Notice that generating a large number of universities can be very time-consuming. To speed up the generation, we advice the usage of a MapReduce program that generates the universities in parallel. This program is available at **https://github.com/jrbn/webpie/blob/master/src/jobs/CreateLUBMDataset.java**.

After a KB has been created (please notice that the triples must be in 'nt' format), launch the command

```
./vlog load -i <path to the nt triples> -o <path to the Trident database>
```

to create a Trident database. Then, modify the **edb.conf** file to point to the newly created Trident database.

Now we can perform materialization using the command

```
./vlog tgchase
```

using the edb.conf file and the "L" Datalog program, as specified above.

### Contributors 

* Efthymia **Tsamoura** (efi.tsamoura@samsung.com), Samsung Research.
* David **Carral** (david.carral@tu-dresden.de), TU Dresden.
* Enrico **Malizia** (enmalizia@gmail.com), University of Bologna.
* Jacopo **Urbani** (jacopo@cs.vu.nl), Vrije Universiteit Amsterdam.


### References

> [1] P. Hu, B. Motik, and I. Horrocks. Modular materialisation of datalog programs. In AAAI, pages 2859–2866, 2019.
