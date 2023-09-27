# /bin/sh

CMP=$PWD/test/compare.py
DATA_PATH=$PWD/data/
CMD=$PWD/src/logis_batch

cd ${DATA_PATH}

gzip -d *.gz

$CMD -B MSigDB.Hallmark.background -X CytoSig.signature -Y MSigDB.Hallmark.mat -out test.output

if [ $? -ne 0 ]         # Test exit status of "logis_batch"
then
  echo "logis_batch failure."
  exit 99
fi

python3 $CMP output test.output
STATUS=$?

rm test.output*
gzip *.*

exit $STATUS
