# ----- USO DE GIT -----
					
### Instalación
En la terminal de UNIX:

		sudo apt-get install git

Luego hay que configurar el nombre y mail del usuario con:

		git config --global <nombre del usuario>
		git config --global <e-mail del usuario>

Despues nos paramos en la carpeta que vamos a trabajar (p.e. TP1) e inicializamos git con:

		git init
	
### Básico

Git tienen tres etapas:

    NORMAL -- INDICE -- COMMIT
	
Todo lo que se haga estará en la etapa NORMAL. Para presentar en el índice (*stage*) y seguir a los archivo(s):

    git add <archivo1> <archivo2>

La primera vez que se agregue el archivo, git seguirá ese archivo. Cada que vez que se modifique un archivo en seguimiento, puede presentarse nuevamente (¡esta vez no lo seguirá porque ya lo hace!). Para ver el estado actual del repositorio utilizar:

    git status

Cuando los archivos estén listos para almacenarse y agregarse a la historia del proyecto (repositorio) escribir en la terminal:

    git commit -m "<descripción del commit>"

Si los archivos a *comitear* ya están en seguimiento, el flag **-a** en commit automáticamente los presenta y commitea en una sola orden.
Para ver la historia del proyecto, donde se enlistan los commit con sus respectivas descripciones, invocar:

    git log

Para ver las diferencias entre los archivos modificados (no *comiteados*) y el último commit en el branch actual:

    git diff
    
Allí detallará con un '-' delante de lo que había antes (en el commit) y con un '+' lo que hay actualmente
Cuando se quiera subir los commit al repositorio remoto utilizar:

    git push -u origin master

Luego pedirá el nombre de usuario y contraseña para finalizar la operación. Es una buena práctica hacer uso del 'rebase' antes de un push para subir el historial más claro.
Para sincronizar los cambios en el repositorio remoto con los propios se usa:

    git pull

Para ilustrar graficamente todo el historial:

    sudo apt-get install gitk
    gitk
    
y para hacer operaciones:

    git gui

### Extra

Para volver atrás con un archivo sin haberlo ingresado al índice (si fue eliminado o si fue modificado, vuelve al estado del último commit):

    git checkout <nombre del archivo>

Si fue agregado pero no commiteado

		git reset HEAD <nombre del archivo>
		git checkout <nombre del archivo>

Si es un directorio:

		git reset HEAD -- <nombre del directorio>

Si fue commiteado (revertir un commit)		[El nombre del commit son los numeros que le siguen cuando hacés git log]

		git revert <nombre del commit>

Si fue agregado al índice por error:

		git reset <nombre del archivo>

-----

Para agregar tags para volver fácilmente a versiones anteriores de commits:

		git tag <nombre del tag> -m "<descripción del tag>"		// Etiqueta el estado actual
		git checkout <nombre del tag>						      // Vuelve todo al estado del tag

Para ver los tags

		git tag

-----

Para ver los branches del repositorio (marca el actual con un '*'):

		git branch

Para crear un branch

		git branch <nombre del branch>

Para cambiar de branch

		git checkout <nombre del branch>

Para eliminar un branch (siempre y cuando no haya conflictos, p.e. con merging)

		git branch -d <nomber del branch>

### MERGE

Merge te permite combinar los cambios de dos branches (¡uno de ellos puede ser remoto!). Merge realiza el comunmente llamado merge de tres versiones entre la última imagen de dos branches, basado en el más reciente ancestro en común de ambos. Como resultado se obtiene una nueva imagen. Es posible hacer un merge entre los cambios de un branch y el actualmente activo mediante el siguiente comando:

		git merge <nombre del branch>

Si un conflicto ocurre, Git marcará el conflicto EN EL MISMO ARCHIVO y deberá resolverse manualmente. Después de resolverlo, se debe agregar el archivo al índice y comitear el cambio. Recién ahí se podrá seguir trabajando. Una herramienta para solucionar los conflictos por merging es:

		git mergetool

------

Para hacer cambios en la historia (cambio en los commits):

		git rebase HEAD~<número de commits que se quieran combinar>

Hacer esto siempre antes de pushear los commits para evitar problemas y errores. Es una buena práctica usar rebase ANTES de pushear para mantener limpio el historial del repositorio remoto.


Para más información visitar:
http://blog.santiagobasulto.com.ar/programacion/2011/11/27/tutorial-de-git-en-espanol.html
https://www.atlassian.com/git/
