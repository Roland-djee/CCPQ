#!/bin/sh

cp ../input/System.inp ../input/System_tmp.inp 
cp ../input/Dynamics.inp ../input/Dynamics_tmp.inp 
cp ../input/Output.inp ../input/Output_tmp.inp 

cp ./features/System.inp ../input/System.inp 
cp ./features/Dynamics.inp ../input/Dynamics.inp 
cp ./features/Output.inp ../input/Output.inp 

./test_read.exe

cp ../input/System_tmp.inp ../input/System.inp 
cp ../input/Dynamics_tmp.inp ../input/Dynamics.inp 
cp ../input/Output_tmp.inp ../input/Output.inp 

rm ../input/System_tmp.inp
rm ../input/Dynamics_tmp.inp
rm ../input/Output_tmp.inp




































































