#!/bin/bash

# usage:
# generate_info_header <compiler> <flags> <header_filename> <macro_prefix>
# example:
#   % generate_info_header g++ -O3 miniFE MINIFE
# this will cause the appropriate info to be put in a 
# header named miniFE_info.hpp and the info will be in macros
# that start with MINIFE.
#
# an example of usage can be seen in miniFE/make_targets
#
if [ $# != 4 ] ; then
    echo "error, need 4 arguments.";
    exit 1;
fi

cxx=`which ${1}`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    cxx="unknown";
fi
echo "CXX: ${cxx}"

cxx_ver=`${1} --version 2>&1`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    cxx_ver=`${1} -V 2>&1`;
    errcode="$?"
    if [ ${errcode} != "0" ] ; then
	cxx_ver="unknown";
    fi
fi

cxx_ver=${cxx_ver// /@}
cxx_version=""
for i in $(echo ${cxx_ver});
do
    if [ "$cxx_version" == "" ]; then
	cxx_version=$i;
    fi
done
cxx_version=${cxx_version//@/ }
echo "Compiler version: ${cxx_version}"

cxxflags=${2}
hostname=`uname -n`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    hostname="unknown";
fi

kern_name=`uname -s`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    kern_name="unknown";
fi

kern_rel=`uname -r`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    kern_rel="unknown";
fi

proc=`uname -p`
errcode="$?"
if [ ${errcode} != "0" ] ; then
    proc="unknown";
fi

header_filename=${3}
macro_prefix=${4}

cat << END_CAT > ${header_filename}
#ifndef ${header_filename/./_}
#define ${header_filename/./_}

#define ${macro_prefix}_HOSTNAME "${hostname}"
#define ${macro_prefix}_KERNEL_NAME "'${kern_name}'"
#define ${macro_prefix}_KERNEL_RELEASE "'${kern_rel}'"
#define ${macro_prefix}_PROCESSOR "'${proc}'"

#define ${macro_prefix}_CXX "'${cxx}'"
#define ${macro_prefix}_CXX_VERSION "'${cxx_version}'"
#define ${macro_prefix}_CXXFLAGS "'${cxxflags}'"

#endif
END_CAT

echo "Generated ${header_filename}"
