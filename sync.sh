#!/bin/bash

usage() {
	echo Should be: sync.sh src dest
	exit 1
}

unis=("МГТУ им. Баумана" МФТИ МГУ)
blocks=(Лекции Семинары .src)

sups=(.addendum .data .fig)

if [ $# -lt 2 ]; then usage; fi

for uni in "${unis[@]}"; do
	for blk in "${blocks[@]}"; do
		if [ ! -d "$1/${uni}/${blk}" ]; then continue; fi
		if [ ! -e "$2/${uni}/${blk}" ]; then mkdir -p "$2/${uni}/${blk}"; fi
		rsync -ahvzu --safe-links --exclude="*/" "$1/${uni}/${blk}/" "$2/${uni}/${blk}/"
	done;
done

for sup in "${sups[@]}"; do
	if [ ! -e "$2/${sup}" ]; then mkdir -p "$2/${sup}"; fi
	rsync -ahvzu --safe-links "$1/${sup}/" "$2/${sup}/"
done
