## Small python library to read bgen format.

This parser is a part of the [HASE](https://github.com/roshchupkin/hase) framework for fast HD GWAS analysis, and provides just basic API for bgen data files reading and manipulation.
Below you can find several examples how you can get data in python format for further analysis.

## Support
* bgen v1.1; 1.2; 1.3
* Layout 1,2

**Fits for UK Biobank data**
Want to convert to more efficient data format? Check  [HASE](https://github.com/roshchupkin/hase)

## Does not support
* Ploidy > 2
* Number of allele > 2
* Phase data

**coming soon ...**

## Installation
1. `git clone  https://github.com/roshchupkin/pybgen.git`
2. Add path to the cloned repository into your python search:
a. `export PYTHONPATH=$PYTHONPATH:{path to pybgen folder}`
b. inside python:
```
>> import sys
>> sys.path.append(path to pybgen folder)
>> import pybgen
 ```

## Requirements

Python library:
1. numpy
2. bitarray
3. zstd (you need this for bgen v1.2). Many thanks!!! to Sergey for simple [python zstd library](https://github.com/sergey-dryabzhinsky/python-zstd)

## Usage

You do not need to have bgen.bgi files. This parser works with pure bgen files and can make its own indices small files.
### Overview
```
>> import pybgen
>> B_test=pybgen.Bgen('example.bgen')
File zise is 665108 bytes
There are 199 variants
There are 500 individuals
Genotype block layout 2

>>  B_test.info()
Name:example.bgen; N samples:500; N probes:199; Compression:zlib; Layout:2

>> B_test.get_indices()
>> B_test.probes_info.keys()[:10]
[u'RSID_2',
 u'RSID_3',
 u'RSID_4',
 u'RSID_5',
 u'RSID_6',
 u'RSID_7',
 u'RSID_8',
 u'RSID_9',
 u'RSID_10',
 u'RSID_11']

>> probe=B_test.read_probe(rsid='RSID_2')
>> probe.info()
Iden: SNPID_2, RSID: RSID_2, CHR: 1, POS: 2000, Alleles: OrderedDict([(1, [u'A']), (2, [u'G'])])

>> probe.get_genotypes(genotypes=True)
>> print probe.prob[:10]
[ 0.          0.          0.02780236  0.00863674  0.01736504  0.04968414
  0.02487179  0.93283081  0.03460688  0.01919559]

>> print probe.genotypes[:10]
[ 0.          0.06424146  0.08441421  0.9825744   0.08840936  0.14108266
  1.07330097  0.05413817  0.10858148  0.12307751]
```

### Make indices file
```
>> B_test.get_indices()
>> B_test.save_indices('/home/username/bgen/')
```

This will save idices files 'example.bgen_ind.npy' to chosen folder.
Nest time you can directly load this info
```
>> B_test=pybgen.Bgen('example.bgen')
>> B_test.load_indices(/home/username/bgen/example.bgen_ind.npy)
```

Actually the operation `get_indices()` does not take a lot of time, but for very intense use of the same bgen files can be quite useful.


## Contacts

If you have any questions/suggestions/comments or problems do not hesitate to contact me!

* gena.roshchupkin@gmail.com
