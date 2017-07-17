For program:
 .yaml files contain parameters; can be created using any text-editor
 .config files contain initial MT configs; created using ./buildmts
 .file files are the output of ./sim

For makefile:
# DEBUG mode
    make sim

# RELEASE mode
	make CFG=release sim
