#!/bin/sh

cp ../Input/System.inp ../Input/System_tmp.inp 
cp ../Input/Dynamics.inp ../Input/Dynamics_tmp.inp 

cp ./features/System.inp ../Input/System.inp 
cp ./features/Dynamics.inp ../Input/Dynamics.inp 

./test_read.exe

cp ../Input/System_tmp.inp ../Input/System.inp 
cp ../Input/Dynamics_tmp.inp ../Input/Dynamics.inp 

rm ../Input/System_tmp.inp
rm ../Input/Dynamics_tmp.inp




































































