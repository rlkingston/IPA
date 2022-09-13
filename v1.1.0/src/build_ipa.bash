
#!/bin/bash
# This is for developers to ensure that the svn-version.h header is updated correctly when building.
# Builds a C++ header containing the SVN version number, the compile date, and the public release version. 

SVN_VER="$(svnversion -n .)"
COMPILE_DATE=`date +%Y/%m/%d`
RELEASE_VER=1.1.0

cat << EOF > svn-version.h
#define xmake_string(x) make_string(x)
#define make_string(x) #x
#define svn_version $SVN_VER
#define compile_date $COMPILE_DATE
#define release_version $RELEASE_VER
EOF

# build ipa from the makefile
make build