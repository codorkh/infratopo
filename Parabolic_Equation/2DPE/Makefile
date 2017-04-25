#=======================================
#
#	   Main Makefile Rules
#
#=======================================

include Make.inc

all: cnpe gtpe gfpe

cnpe:
	(cd SRC/CNPE; make)

gtpe:
	(cd SRC/GTPE; make)

gfpe:
	(cd SRC/GFPE; make)

clean:
	(cd SRC/CNPE; make clean)
	(cd SRC/GTPE; make clean)
	(cd SRC/GFPE; make clean)
