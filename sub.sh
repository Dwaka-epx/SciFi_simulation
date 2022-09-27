#!/bin/bash 

# https://yutarine.blogspot.com/2018/02/shell-command-bc02.html
# [コマンド] bc(数値計算ソフト)の結果出力で、小数点前の0(ゼロ)を省略させないようにする方法

processes=(
build_elec
build_kaon
build_muon
build_pion
build_prot
)

# loop over processes
for process in "${processes[@]}"; do
	cd ${process}
	#rm -rf CMake*; ./compile.sh
   make -j5
   cp -r ../run.sh ./
	./run.sh
   cd ..	
done

