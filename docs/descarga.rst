descarga module
===============
Para obtener el ADN que se utilizará como referencia para el diseño de las sondas se va a descargar directamente el archivo genético desde la base de datos ‘Nucleotide’ de NCBI (National Center for Biotechnology Information) del NIH (National Institute of Health, EEUU). Esta base de datos contiene millones de registros de secuencias de nucleótidos con mucha información adicional, como la especie, cromosoma, nombre del gen, exones, la fuente de los datos, etc.

Para llevar a cabo la descarga de los archivos genéticos, se implementó un módulo que permite realizar descargas desde esta base de datos a través de la librería Entrez y almacenar los archivos en una carpeta. De esta manera se puede automatizar la obtención de las secuencias de referencia.

Por otra parte se implementó la función que recibe un archivo almacenado localmente y lo retorna como SeqRecord, además de crear una copia del archivo y almacenarla. Esta función sólo podrá recibir archivos de secuencia anotados como GenBank y GFF3.

El producto de estas funciones implementadas es la secuencia en formato SeqRecord, una estructura de datos propia de la librería BioPython que permite almacenar grandes secuencias con su respectiva notación, como por ejemplo las posiciones de cada uno de los exones, las transcripciones, secciones no codificantes, entre otros. Además contiene información extra de la secuencia tal como la especie, el cromosoma, la fuente de donde salieron los datos, etc.

.. automodule:: descarga
   :members:
   :undoc-members:
   :show-inheritance:
