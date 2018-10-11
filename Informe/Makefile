# Makefile para informes en LaTex utilizando además los programas Octave y XCircuit

# Se supone que el Makefile se encuentra en una carpeta LaTeX con los .tex, mientras
# que la carpeta LaTeX se encuentra en una carpeta mayor que también contiene los otros
# programas con sus respectivas carpetas. El informe se imprime en la carpeta ppal del 
# proyecto

## .PHONY se utiliza para destacar que las reglas con targets 'phoneados' no generan 
## archivos con esos nombres. Por lo tanto, no busca la existencia de ellos y opera solo
## si es llamado y siempre. También se usa para evitar ambiguedades (como con el comando view)


NAME_PROJECT = TP1 
TEX = pdflatex
SEC = ./Secciones/*
FIG = ./Figuras/*
OUTPUT_DIR = ../
FIG_DIR = ./Figuras
GRAPHICS_DIR = ../Octave
CIRCUITS_DIR = ../XCircuit


all: main.pdf rename_copy clean
release: main.pdf cls


main.pdf: main.tex config.tex preamble.tex $(SEC) $(FIG)
	$(TEX) -halt-on-error main.tex
	$(TEX) main.tex
# -halt-on-error: Cuando hay un error en la compilación corta sin preguntar qué hacer.


.PHONY: rename_copy
rename_copy:
	cp main.pdf $(OUTPUT_DIR)$(NAME_PROJECT).pdf
# Se copia en vez de moverse, porque Make ve que falta lo generado por main y compila de nuevo
# // Podría hacerse en el make con "-jobname <name>" y "-output-directory <dir>" pero habría
# // que borrar los .aux .log y .out en ../


.PHONY: view
view:
	evince $(OUTPUT_DIR)$(NAME_PROJECT).pdf &


.PHONY: clean
clean:
	-rm -f *.out *.aux *.log
	-rm -f *.lof *.lot *.toc *.lol
	-rm -f *~
	-rm -f ./Secciones/*~

.PHONY: cls
cls: clean
	-rm -f $(FIG_DIR)/*eps-converted-to.pdf
	-rm -f ./main.pdf
#	make -C $(GRAPHICS_DIR) clean
#	make -C $(CIRCUITS_DIR) clean
# El -rm elimina el archivo (rm) pero sin quejarse si no puede borrarlo (-).
# El -f no pone texto si de que se pudo o no pudo.

## El clean borra en cada compilación mientras que el cls se encarga de borrar
## al finalizar el proyecto todo lo que puede ser generado con los scripts.



.PHONY: graphics
graphics:
	make -C $(GRAPHICS_DIR) move
	make all
# Si pusiese all en vez de print, quedo dentro del Octave y puedo trabajar allí


.PHONY: circuits
circuits:
	make -C $(CIRCUITS_DIR) all
	make all


# Al finalizar el proyecto, se eliminan todos los archivos autogenerados para mejor almacenamiento

