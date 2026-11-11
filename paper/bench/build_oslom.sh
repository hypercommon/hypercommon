set -e
mkdir -p ~/tools && cd ~/tools
curl -sL -o OSLOM2.tar.gz http://oslom.org/code/OSLOM2.tar.gz
ls -la OSLOM2.tar.gz
tar xzf OSLOM2.tar.gz
cd OSLOM2
CXX=$(command -v g++ || ls ~/tools/mamba/envs/gt/bin/*-g++ 2>/dev/null | head -1)
echo "compiler $CXX"
[ "$(uname)" = Darwin ] || sed -i "s#^g++#$CXX#; s# g++ # $CXX #g" compile_all.sh || true
sh compile_all.sh > compile.log 2>&1 || (tail -20 compile.log; exit 1)
ls -la ~/tools/OSLOM2/oslom_undir
echo "OSLOM BUILD OK"
