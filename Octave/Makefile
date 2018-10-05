# Makefile para scripts de Octave

SCRIPT_NAME = graph.m
PROGRAM = octave
DEBUG = -d
OUTPUT_DIR = ./
FIG_DIR = ../Informe/Figuras
#SPICE_DIR = ../LTSpice
#LAB_DIR = ../Laboratorio



all: 
	$(PROGRAM) --persist $(SCRIPT_NAME) -p 0
	-rm -f octave-workspace
	-rm -f *~

debug: 
	$(PROGRAM) $(DEBUG) $(SCRIPT_NAME)


print:
	$(PROGRAM) -q $(SCRIPT_NAME) -p 1
#	mv $(OUTPUT_DIR)*.eps $(FIG_DIR)
#	mv $(OUTPUT_DIR)*.tex $(FIG_DIR)
	-rm -f octave-workspace
	-rm -f *~

spice:
	make -C $(SPICE_DIR) all
# Copia los .txt con los datos

lab:
	make -C $(LAB_DIR) all
# Copia los .txt con las mediciones

.PHONY: move
move:
	mv $(OUTPUT_DIR)*.eps $(FIG_DIR)
	mv $(OUTPUT_DIR)*.tex $(FIG_DIR)

.PHONY: clean
clean:
	-rm -f octave-workspace
	-rm -f *~
	make -C $(SPICE_DIR) clean
	make -C $(LAB_DIR) clean

ppython:
	py graficar.py
