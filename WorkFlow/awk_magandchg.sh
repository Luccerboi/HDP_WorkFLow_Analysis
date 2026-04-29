#!/bin/bash

mB1=$(awk '/magnetization \(x\)/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==1 {print $NF}')
mB2=$(awk '/magnetization \(x\)/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==2 {print $NF}')
mA=$(awk '/magnetization \(x\)/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==3 {print $NF}')
mX=$(awk '/magnetization \(x\)/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==5 {print $NF}')
mtot=$(awk '/magnetization \(x\)/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1=="tot" {print $NF}')

chgB1=$(awk '/total charge$/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==1 {print $NF}')
chgB2=$(awk '/total charge$/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==2 {print $NF}')
chgA=$(awk '/total charge$/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==3 {print $NF}')
chgX=$(awk '/total charge$/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1==5 {print $NF}')
chgtot=$(awk '/total charge$/,/^tot/ {print$0}' OUTCAR | awk '!a[$0]++' | awk '$1=="tot" {print $NF}')
ZTOT=$(grep ZVAL POTCAR | awk 'BEGIN{ZTOT=0}; {print $6}')

echo $mB1 $mB2 $mA $mX  $mtot

echo $chgB1 $chgB2 $chgA $chgX $chgtot 

echo $ZTOT