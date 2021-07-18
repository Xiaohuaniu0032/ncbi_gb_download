cd HBV
split -l 10000 ../../HBV.list.seq -d -a 4 list_

cd ../HCV
split -l 10000 ../../HCV.list.seq -d -a 4 list_
