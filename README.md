## S<sup>3</sup> v1.0 ##

***S3*** (Sequence Similarity Search) is a program for searching DNA or peptide sequences in protein database.


### Installation ###
(1) Download the source code from Github

    git clone git@github.com:acgtun/s3.git

(2) Build and Install
    
    cd s3
    make all
    make install


### Indexing Genome ###
    
    makedb -d <protein database> -o <index file>

### Searching ###

    s3 -i <index file> -q <read files> -o <output file> [options]



### Searching Options ###


| Option | Long Tag | Type | Default | Brief Description |
| :-------------: |:-------------:|:-----:|:-----:| :-----|
| -i      | -index | String | NULL | index file created by ***makedb*** command ( .dbindex) |
| -q      | -query | String | NULL | query file (.fastq or .fq) |

To see the list of options, use "-?" or "-help".

    
### Contacts ###

***Haifeng Chen***   *haifengc@usc.edu*

***Ting Chen***   *tingchen@usc.edu*


### Copyright ###

    Copyright (C) 2016 University of Southern California
                   
    Authors: Haifeng Chen and Ting Chen

    S3 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    S3 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with S3.  If not, see <http://www.gnu.org/licenses/>.
