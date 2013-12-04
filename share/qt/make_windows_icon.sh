#!/bin/bash
# create multiresolution windows icon
ICON_SRC=../../src/qt/res/icons/datacoin.png
ICON_DST=../../src/qt/res/icons/datacoin.ico
convert ${ICON_SRC} -resize 16x16 datacoin-16.png
convert ${ICON_SRC} -resize 32x32 datacoin-32.png
convert ${ICON_SRC} -resize 48x48 datacoin-48.png
convert ${ICON_SRC} -resize 64x64 datacoin-64.png
convert datacoin-32.png ${ICON_SRC} datacoin-64.png datacoin-48.png datacoin-32.png datacoin-16.png ${ICON_DST}

