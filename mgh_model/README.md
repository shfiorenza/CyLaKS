For program:
 .yaml files contain parameters; can be created using any text-editor
 .file files are the output of ./sim

For makefile:
# DEBUG mode
    make sim

# RELEASE mode (optimized compiler flags; MUCH faster runtime)
	make CFG=release sim
