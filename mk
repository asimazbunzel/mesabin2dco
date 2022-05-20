#!/bin/bash

function check_okay {
	if [ $? -ne 0 ]
	then
		echo
		echo "FAILED"
		echo
		exit 1
	fi
}

function mkdir_mods {
   if [ ! -r lib ]
   then
      mkdir lib
   fi

   if [ ! -r include ]
   then
      mkdir include
   fi
}

function do_one {
   cd src/$1
   check_okay
   echo
   pwd
   echo "building $1 package"
   echo
   ./mk
   check_okay
   ./export
   check_okay
   echo "mesa/$1 has been built and exported"
   echo
   cd ../..
}

# make dirs
mkdir_mods; check_okay

# install modules
do_one ce
do_one core_collapse

# link everything with the source code
cd make; make; check_okay
